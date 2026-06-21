from lifton import intervals, lifton_utils, annotation, extract_sequence, stats, logger, run_liftoff, run_miniprot, run_evaluation, gff3_validator, __version__
from intervaltree import Interval
import argparse
from pyfaidx import Fasta
import os, sys
import time
import concurrent.futures


def _describe_annotation_source(x):
    """Render a Liftoff/miniprot annotation DB source for the
    "Creating X annotation database" log lines.

    Normally a path string. Under ``--inmemory-liftoff`` and
    ``--stream``, however, the value is the full GFF3 as an in-memory
    ``bytes`` blob. f-string interpolation of that blob produces a
    single multi-megabyte stderr line (Phase 16 follow-up: the bee
    benchmark previously emitted a 33 MB single-line stderr from
    this site). Bytes values are rendered as a brief length summary
    so log files stay grep-able and small; path values pass through
    unchanged.
    """
    if isinstance(x, (bytes, bytearray)):
        return f"<in-memory bytes, {len(x):,} bytes>"
    return x

def args_gffutils(parser):
    gffutils_grp = parser.add_argument_group('* gffutils parameters')
    gffutils_grp.add_argument(
        '--merge-strategy', default='create_unique', choices=['create_unique', 'merge', 'error', 'warning', 'replace'],
        help='Strategy for merging features when building the database. Default is "create_unique".'
    )
    gffutils_grp.add_argument(
        '--id-spec', default=None,
        help='Attribute to use as feature ID. Default is "ID".'
    )
    gffutils_grp.add_argument(
        '--force', default=False, action='store_true', help='Overwrite existing database.'
    )
    gffutils_grp.add_argument(
        '--verbose', default=False, action='store_true', help='Enable verbose output.'
    )
    return gffutils_grp


def args_outgrp(parser):
    outgrp = parser.add_argument_group('* Output settings')
    outgrp.add_argument(
        '-o', '--output', default='lifton.gff3', metavar='FILE',
        help='write output to FILE in same format as input; by default, output is written to "lifton.gff3"'
    )
    outgrp.add_argument(
        '-u', default='unmapped_features.txt', metavar='FILE',
        help='write unmapped features to FILE; default is "unmapped_features.txt"',
    )
    outgrp.add_argument(
        '-exclude_partial', action='store_true',
        help='write partial mappings below -s and -a threshold to unmapped_features.txt; if true '
             'partial/low sequence identity mappings will be included in the gff file with '
             'partial_mapping=True, low_identity=True in comments'
    )
    return outgrp


def args_aligngrp(parser):
    aligngrp = parser.add_argument_group('Alignments')
    aligngrp.add_argument('-mm2_options', metavar='=STR', type=str, default='-a --end-bonus '
                                                                            '5 --eqx -N 50 '
                                                                            '-p 0.5',
                          help='space delimited minimap2 parameters. By default ="-a --end-bonus 5 --eqx -N 50 -p 0.5"')
    aligngrp.add_argument('-mp_options', metavar='=STR', type=str, default='',
                          help='space delimited miniprot parameters.')
    aligngrp.add_argument(
        '-a', default=0.5, metavar='A', type=float,
        help='designate a feature mapped only if it aligns with coverage ≥A; by default A=0.5',
    )
    aligngrp.add_argument(
        '-s', default=0.5, metavar='S', type=float,
        help='designate a feature mapped only if its child features (usually exons/CDS) align '
             'with sequence identity ≥S; by default S=0.5'
    )
    aligngrp.add_argument(
        '-min_miniprot', default=0.9, metavar='MIN_MINIPROT', type=float,
        help='The minimum length ratio of a protein-coding transcript to the longest protein-coding transcript within a gene locus, as identified exclusively by miniprot in the target genome, is set by default to MIN_MINIPROT=0.9.'
    )
    aligngrp.add_argument(
        '-max_miniprot', default=1.5, metavar='MAX_MINIPROT', type=float,
        help='The maximum length ratio of a protein-coding transcript to the longest protein-coding transcript within a gene locus, as identified exclusively by miniprot in the target genome, is set by default to MIN_MINIPROT=1.5.'
    )
    aligngrp.add_argument(
        '-d', metavar='D', default=2.0, type=float,
        help='distance scaling factor; alignment nodes separated by more than a factor of D in '
             'the target genome will not be connected in the graph; by default D=2.0'
    )
    aligngrp.add_argument(
        '-flank', default=0, metavar='F', type=float, help="amount of flanking sequence to align as a "
                                                           "fraction [0.0-1.0] of gene length. This can improve gene "
                                                           "alignment where gene structure  differs between "
                                                           "target and "
                                                           "reference; by default F=0.0")
    return aligngrp


def args_optional(parser):
    parser.add_argument('-V', '--version', help='show program version', action='version', version=__version__)
    parser.add_argument('-D', '--debug', action='store_true', help='Run debug mode')
    parser.add_argument(
        '-t', '--threads', default=1, type=int, metavar='THREADS', help='use t parallel processes to accelerate alignment; by default p=1'
    )
    parser.add_argument('-m', help='Minimap2 path', metavar='PATH')
    parser.add_argument('-f', '--features', metavar='TYPES', help='list of feature types to lift over')
    parser.add_argument(
        '-infer-genes', required=False, action='store_true',
        help='use if annotation file only includes transcripts, exon/CDS features. '
             'Automatically enabled for GTF files.'
    )
    parser.add_argument(
        '-infer_transcripts', action='store_true', required=False,
        help='use if annotation file only includes exon/CDS features and does not include transcripts/mRNA. '
             'Automatically enabled for GTF files.'
    )
    parser.add_argument(
        '-chroms', metavar='TXT', help='comma seperated file with corresponding chromosomes in '
                                       'the reference,target sequences',
    )
    parser.add_argument(
        '-unplaced', metavar='TXT',
        help='text file with name(s) of unplaced sequences to map genes from after genes from '
             'chromosomes in chroms.txt are mapped; default is "unplaced_seq_names.txt"',
    )
    parser.add_argument('-copies', action='store_true', help='look for extra gene copies in the target genome')
    parser.add_argument(
        '-sc', default=1.0, metavar='SC', type=float,
        help='with -copies, minimum sequence identity in exons/CDS for which a gene is considered '
             'a copy; must be greater than -s; default is 1.0',
    )
    parser.add_argument('-overlap', default=0.1, metavar='O', help="maximum fraction [0.0-1.0] of overlap allowed by 2 "
                                                                   "features; by default O=0.1", type=float)
    parser.add_argument('-mismatch', default=2, metavar='M', help="mismatch penalty in exons when finding best "
                                                                  "mapping; by default M=2", type=int)
    parser.add_argument('-gap_open', default=2, metavar='GO', help="gap open penalty in exons when finding best "
                                                                   "mapping; by default GO=2", type=int)
    parser.add_argument('-gap_extend', default=1, metavar='GE', help="gap extend penalty in exons when finding best "
                                                                     "mapping; by default GE=1", type=int)
    parser.add_argument('-subcommand', required=False,  help=argparse.SUPPRESS)
    parser.add_argument('-polish', required=False, action='store_true', default = False)
    parser.add_argument('-cds', required=False, action="store_true", default=True, help="annotate status of each CDS "
                                                                                        "(partial, missing start, "
                                                                                        "missing stop, inframe stop "
                                                                                        "codon)")
    parser.add_argument(
        '-time', '--measure_time', required=False, action='store_true',
        help='Enable time measurement for each step'
    )
    parser.add_argument(
        '--validate-output', required=False, action='store_true', default=False,
        help='Validate the generated GFF3 output file for format correctness and '
             'feature hierarchy after writing. Prints a detailed report to stderr.'
    )
    parser.add_argument(
        '--validate-verbose', required=False, action='store_true', default=False,
        help='When --validate-output is set, also print warnings (not just errors)'
    )
    parser.add_argument(
        '--strict-gff', dest='strict_gff', action='store_true', default=False,
        help='Run the NCBI GFF3 input-side validator on the reference '
             'annotation and exit non-zero on any spec violation '
             '(missing ##gff-version 3, start>end, negative coords, '
             'unencoded reserved chars, dangling Parent, etc.).'
    )
    parser.add_argument(
        '--stream', dest='stream', action='store_true', default=False,
        help='Phase 7 streaming-adapter fast path: pipe miniprot stdout '
             'directly into an in-memory gffbase FeatureDB instead of '
             'writing miniprot.gff3 to disk. Eliminates the SQLite '
             're-ingest of miniprot output. Output GFF3 is byte-identical '
             'to the default path; this flag changes I/O, not algorithms.'
    )
    parser.add_argument(
        '--inmemory-liftoff', dest='inmemory_liftoff', action='store_true',
        default=False,
        help='Phase 8 in-memory Liftoff fast path: serialise Liftoff\'s '
             'lifted_feature_list to bytes inside the parent process and '
             'feed it straight to gffbase, skipping the liftoff.gff3 disk '
             'write and SQLite re-ingest. Output GFF3 is byte-identical '
             'to the default path; this flag changes I/O, not algorithms.'
    )
    parser.add_argument(
        '--locus-pipeline', dest='locus_pipeline', action='store_true',
        default=False,
        help='Locus-major fan-out: dispatch Step 7 (per-Liftoff-gene '
             'processing) through a ThreadPoolExecutor sized by --threads. '
             'Output is emitted in submission order so --threads N is '
             'byte-identical to --threads 1; this flag changes scheduling, '
             'not algorithms. As of Iteration 8 this works on the DEFAULT '
             '(gffutils) backend WITHOUT --native — per-locus work runs '
             'against materialised proxy DBs, so any backend is thread-safe '
             '(set LIFTON_PARALLEL_BLOCK_GFFUTILS=1 to opt back out to '
             'serial-on-gffutils). Note: combined with the default '
             'concurrent Step 4, peak busy cores can reach ~N+1.'
    )
    parser.add_argument(
        '--native', dest='native', action='store_true', default=False,
        help='Phase 10 native bindings: route miniprot through the '
             'pyminiprot-shaped facade and unlock in-process Step-7 '
             'threading (Phase 17b). As of Iteration 7, --native does NOT '
             'route Liftoff\'s minimap2 alignment through `mappy` by '
             'default — that in-process path is slower and slightly less '
             'accurate than the subprocess minimap2 path (mappy can\'t take '
             '--end-bonus/-p), so Liftoff alignment stays on the proven '
             'subprocess path. Set LIFTON_NATIVE_LIFTOFF_ALIGN=1 to opt the '
             'mappy Liftoff path back in (e.g. no minimap2 binary on PATH). '
             'Falls back gracefully when mappy is not installed.'
    )
    parser.add_argument(
        '--serial-aligners', dest='serial_aligners', action='store_true',
        default=False,
        help='Restore the pre-Iteration-6 SEQUENTIAL Step 4: run Liftoff '
             '(DNA) then miniprot (protein) one after the other instead of '
             'the default (concurrent) overlap. The default now dispatches '
             'miniprot (an independent subprocess) to a background thread '
             'while Liftoff runs on the main thread (Liftoff reads the '
             'main-thread-bound SQLite reference DB, so it cannot move off '
             'it), collapsing Step-4 wall from t_liftoff + t_miniprot to '
             'max(t_liftoff, t_miniprot). The concurrent default is '
             'byte-identical to this serial path (only miniprot\'s timing '
             'moves); use --serial-aligners on core-constrained machines '
             '(concurrent peak is ~N+1 cores with --threads N) or to keep '
             'the two tools\' console logs from interleaving.'
    )
    parser.add_argument(
        '--parallel-aligners', dest='parallel_aligners', action='store_true',
        default=False,
        help='No-op alias (kept for backward compatibility). The Step-4 '
             'Liftoff/miniprot overlap that this flag used to gate is now '
             'the DEFAULT (Iteration 6 promotion), so --parallel-aligners '
             'has no effect; pass --serial-aligners to opt out.'
    )
    parser.add_argument(
        '--optimize', dest='optimize', action='store_true', default=False,
        help='No-op alias (kept for backward compatibility). The best-of-'
             'outcome verified Liftoff/miniprot merge that this flag used to '
             'gate is now the DEFAULT, so --optimize has no effect. Reserved '
             'for future opt-in v2 accuracy/speed work.'
    )
    parser.add_argument(
        '--legacy-merge', dest='legacy_merge', action='store_true', default=False,
        help='Restore the pre-promotion Liftoff/miniprot merge: apply the '
             'protein-maximization chained CDS UNCONDITIONALLY (no best-of-'
             'outcome verification). This is the published-manuscript merge '
             'and can silently frameshift downstream CDS on divergent inputs; '
             'use only for reproducing legacy output. The default path now '
             'runs the verified best-of-outcome merge instead.'
    )
    parser.add_argument(
        '--full-dp-align', dest='full_dp_align', action='store_true', default=False,
        help='Restore the exact pre-Iteration-3 giant-only alignment: full '
             'Needleman-Wunsch DP for every non-giant gene (gate 8000 aa / '
             '25000 nt); giants still memory-bounded-windowed so it cannot OOM. '
             'The default is now "band everything" (anchor-windowed above '
             '~2500 aa / 8000 nt: 1.4-2.6x faster, identity-exact on '
             'same-species lifts). Use this for manuscript reproduction or '
             'maximal accuracy on divergent inputs. For PURE full DP including '
             'giants, set LIFTON_ALIGN_WINDOW_AA/NT to a huge value.'
    )
    parser.add_argument(
        '--fast-align', dest='fast_align', action='store_true', default=False,
        help='No-op alias (kept for backward compatibility). Band-everything '
             'alignment that this flag used to gate is now the DEFAULT, so '
             '--fast-align has no effect; use --full-dp-align to opt OUT.'
    )
    parser.add_argument(
        '--gene-only', dest='gene_only', action='store_true', default=False,
        help='Restore the pre-Iteration-12 GENE-ONLY lift: process only the '
             '`gene` hierarchy, dropping every other gene-like top-level '
             'parent type (pseudogenes, ncRNA_genes, structured mobile '
             'elements). The default now lifts ALL auto-detected gene-like '
             'types. Use --gene-only for manuscript reproduction of the old '
             'default, or when you want strictly the `gene` partition. Ignored '
             'if -f/--features is given (your explicit file always wins).'
    )
    parser.add_argument(
        '--lift-gene-like', dest='lift_gene_like', action='store_true', default=False,
        help='No-op alias (kept for backward compatibility). Lifting all '
             'auto-detected gene-like parent types (pseudogenes, ncRNA_genes, '
             'structured mobile elements) beyond `gene` is now the DEFAULT '
             '(Iteration 12 promotion), so --lift-gene-like has no effect; pass '
             '--gene-only to opt OUT and restore the gene-only lift.'
    )
    parser.add_argument(
        '--miniprot-rescue', dest='miniprot_rescue', action='store_true', default=False,
        help='[EXPERIMENTAL] Regime-gated miniprot-only rescue (OFF by default). '
             'For a reference coding gene the DNA lift MISSED ENTIRELY (its '
             'miniprot mRNA overlaps no lifted gene locus), emit the miniprot '
             'model even when its length ratio falls outside the default '
             '-min_miniprot/-max_miniprot band -- gated instead by a '
             'protein-identity floor (LIFTON_MINIPROT_RESCUE_MIN_ID, default '
             '0.5) within a wider sanity band (LIFTON_MINIPROT_RESCUE_LEN=lo,hi, '
             'default 0.5,2.0). Recovers genuinely-missing genes at large '
             'evolutionary distance (net +recall on the distant/very-distant '
             'tier; see benchmarks/compare/miniprot_rescue_ab.md). CAVEAT: ON is '
             'NOT a strict superset of OFF -- an emitted rescue is added to the '
             'Step-8 suppression tree, so two overlapping missing-gene '
             'candidates can SWAP which is emitted and a multi-hit ref '
             'transcript can get a redundant model. A clean separate-pass '
             'version (lifton2-style dedup) is the planned default-ready form. '
             'Env LIFTON_MINIPROT_RESCUE=1 also enables it.'
    )


def parse_args(arglist):
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    parser.add_argument('target', help='target fasta genome to lift genes to')
    parser.add_argument('reference', help='reference fasta genome to lift genes from')
    parser.add_argument('-E', '--evaluation', help='Run LiftOn in evaluation mode', action='store_true', default = False)
    parser.add_argument('-EL', '--evaluation-liftoff-chm13', help='Run LiftOn in evaluation mode', action='store_true', default = False)
    parser.add_argument('-c', '--write_chains', help='Write chaining files', action='store_true', default = True)
    parser.add_argument('--no-orf-search', help='Do not perform open reading frame search', action='store_true', default = False)
    parser_outgrp = args_outgrp(parser)
    parser_aligngrp = args_aligngrp(parser)
    args_optional(parser)
    referencegrp = parser.add_argument_group('* Required input (Reference annotation)')    
    referencegrp.add_argument(
        '-g', '--reference-annotation', metavar='GFF',  required=True,
        help='the reference annotation file to lift over in GFF3 or GTF format (or) '
                'name of feature database. GTF files are automatically detected and converted to GFF3 for better compatibility. '
                'For best results with GTF files, ensure agat or gffread is installed.'
    )
    ###################################
    # START for the LiftOn params
    ###################################
    referenceseqgrp = parser.add_argument_group('* Optional input (Reference sequences)')
    referenceseqgrp.add_argument(
        '-P', '--proteins', metavar='FASTA', required=False, default=None,
        help='the reference protein sequences. If not provided, the protein sequences will be extracted from the reference annotation. The ID of the protein sequences must match the ID of the transcript sequences.'
    )
    referenceseqgrp.add_argument(
        '-T', '--transcripts', metavar='FASTA', required=False, default=None,
        help='the reference transcript sequences. If not provided, the transcript sequences will be extracted from the reference annotation.'
    )
    liftoffrefrgrp = parser.add_argument_group('* Optional input (Liftoff annotation)')
    liftoffrefrgrp.add_argument(
        '-L', '--liftoff', metavar='gff', default=None,
        help='the annotation generated by Liftoff (or) '
                'name of Liftoff gffutils database; if not specified, '
                'a Liftoff database will be built automatically'
    )
    miniprotrefrgrp = parser.add_argument_group('* Optional input (miniprot annotation)')
    miniprotrefrgrp.add_argument(
        '-M', '--miniprot', metavar='gff', default=None,
        help='the annotation generated by miniprot (or) '
                'name of miniprot gffutils database; if not specified, '
                'a miniprot database will be built automatically'
    )
    parser_gffutils_grp = args_gffutils(parser)
    parser.add_argument('-ad', '--annotation-database', metavar='SOURCE', help='The source of the reference annotation (RefSeq / GENCODE / others).', default = "RefSeq")
    parser.add_argument('--no-auto-convert-gtf', action='store_true', default=False,
                        help='Disable automatic GTF to GFF3 conversion. By default, LiftOn will attempt to convert GTF files to GFF3 for better compatibility.')
    ###################################
    # END for the LiftOn params
    ###################################
    parser._positionals.title = '* Required input (sequences)'
    parser._optionals.title = '* Miscellaneous settings'
    parser._action_groups = [parser._positionals, referencegrp, referenceseqgrp, 
                             liftoffrefrgrp, miniprotrefrgrp, parser_gffutils_grp, 
                             parser_outgrp, parser._optionals, parser_aligngrp]
    args = parser.parse_args(arglist)
    if '-a' not in args.mm2_options:
        args.mm2_options += ' -a'
    if '--eqx' not in args.mm2_options:
        args.mm2_options += ' --eqx'
    if '-N' not in args.mm2_options:
        args.mm2_options += " -N 50"
    if '-p' not in args.mm2_options:
        args.mm2_options += " -p 0.5"
    if '--end-bonus' not in args.mm2_options:
        args.mm2_options += " --end-bonus 5"
    if (float(args.s) > float(args.sc)):
        parser.error("-sc must be greater than or equal to -s")
    if (args.chroms is None and args.unplaced is not None):
        parser.error("-unplaced must be used with -chroms")    
    return args
    

def resolve_miniprot_rescue_args(args):
    """Iteration 22: resolve the regime-gated miniprot-only-rescue flag + its
    tunables ONCE, honouring the env escape hatches, and stash them on `args` so
    run_miniprot.process_miniprot (Step 8) can read them. Flag OFF => Step 8 is
    byte-identical to the pre-Iteration-22 path. Pure/idempotent (no I/O beyond
    env reads) so it is unit-testable.

    - LIFTON_MINIPROT_RESCUE=1|true|yes   -> force the flag ON
    - LIFTON_MINIPROT_RESCUE_MIN_ID=<f>   -> protein-identity floor (default 0.5)
    - LIFTON_MINIPROT_RESCUE_LEN=lo,hi    -> wider sanity band (default 0.5,2.0)
    """
    rescue = getattr(args, "miniprot_rescue", False)
    if os.environ.get("LIFTON_MINIPROT_RESCUE", "").lower() in ("1", "true", "yes"):
        rescue = True
    args.miniprot_rescue = rescue
    try:
        args.miniprot_rescue_min_id = float(os.environ.get("LIFTON_MINIPROT_RESCUE_MIN_ID", "0.5"))
    except (ValueError, TypeError):
        args.miniprot_rescue_min_id = 0.5
    try:
        lo, hi = (float(x) for x in os.environ.get("LIFTON_MINIPROT_RESCUE_LEN", "0.5,2.0").split(","))
        args.miniprot_rescue_len = (lo, hi)
    except (ValueError, TypeError):
        args.miniprot_rescue_len = (0.5, 2.0)
    return args


def run_all_lifton_steps(args):
    t1 = time.process_time()
    # Iteration-3 "band everything" alignment is the DEFAULT (set at align-module
    # import). --full-dp-align restores the exact pre-Iteration-3 giant-only path
    # for manuscript reproduction / paranoid accuracy. Lazy import keeps `align`
    # (which pulls in the lifton_class↔lifton_utils cycle) off the module top.
    if getattr(args, "full_dp_align", False):
        from lifton import align as _align
        _align.configure_alignment(band=False)
    resolve_miniprot_rescue_args(args)
    ################################
    # Step 0: Reading target & reference genomes
    ################################
    tgt_genome = args.target
    ref_genome = args.reference
    if args.output == "stdout": 
        outdir = "."
    else:
        outdir = os.path.dirname(args.output)
        outdir = outdir if outdir != "" else "."
        os.makedirs(outdir, exist_ok=True)
    lifton_outdir = f"{outdir}/lifton_output/"
    args.directory = "intermediate_files/"
    intermediate_dir = f"{outdir}/lifton_output/{args.directory}"
    os.makedirs(intermediate_dir, exist_ok=True)
    stats_dir= f"{outdir}/lifton_output/stats/"
    os.makedirs(stats_dir, exist_ok=True)
    args.directory = intermediate_dir
    logger.log(">> Reading target genome ...", debug=True)
    if not os.path.exists(tgt_genome):
        logger.log_error(f"Target genome file not found: {tgt_genome}")
        sys.exit(1)
    try:
        tgt_fai = Fasta(tgt_genome)
    except Exception as e:
        logger.log_error(f"Failed to read/index target genome '{tgt_genome}': {e}")
        sys.exit(1)
        
    logger.log(">> Reading reference genome ...", debug=True)
    if not os.path.exists(ref_genome):
        logger.log_error(f"Reference genome file not found: {ref_genome}")
        sys.exit(1)
    try:
        ref_fai = Fasta(ref_genome)
    except Exception as e:
        logger.log_error(f"Failed to read/index reference genome '{ref_genome}': {e}")
        sys.exit(1)

    t2 = time.process_time()
    ################################
    # Phase 5 bug fix #6: NCBI GFF3 validation gate
    ################################
    from lifton.io.gff3_validator import GFF3Validator
    target_seqids = set(tgt_fai.keys()) | set(ref_fai.keys())
    findings = GFF3Validator(
        target_seqids=target_seqids,
        strict=getattr(args, "strict_gff", False),
    ).validate(args.reference_annotation)
    # Phase 16 Tier 4: real-world NCBI/RefSeq inputs trigger hundreds of
    # thousands of `unencoded_reserved_char` findings on Dbxref values
    # (DBTAG:ID is technically reserved-char-bearing). The previous
    # unconditional per-finding stderr dump produced 100+ MB stderr
    # logs that buried real errors. Strict mode keeps per-row stderr
    # output (users opted in); the default path now writes findings to
    # a side-car file under stats/ and prints one summary line.
    _strict_gff = getattr(args, "strict_gff", False)
    if _strict_gff:
        for f in findings:
            logger.log(str(f), debug=True)
    elif findings:
        findings_path = os.path.join(stats_dir, "gff3_input_validation.txt")
        try:
            with open(findings_path, "w") as fw:
                for f in findings:
                    fw.write(str(f) + "\n")
            n_err = sum(1 for f in findings if f.severity == "error")
            n_warn = len(findings) - n_err
            logger.log_info(
                f">> GFF3 input validator: {len(findings)} finding(s) "
                f"({n_err} error, {n_warn} warning) written to "
                f"{findings_path}; pass --strict-gff to also dump "
                f"per-row to stderr."
            )
        except OSError as exc:
            logger.log_warning(
                f"GFF3 input validator: could not write findings to "
                f"{findings_path}: {exc}; falling back to stderr dump."
            )
            for f in findings:
                logger.log(str(f), debug=True)
    if _strict_gff and any(f.severity == "error" for f in findings):
        sys.exit(2)
    ################################
    # Step 1: Building database from the reference annotation
    ################################
    logger.log("\n>> Creating reference annotation database : ", args.reference_annotation, debug=True)
    auto_convert_gtf = not args.no_auto_convert_gtf
    ref_db = annotation.Annotation(args.reference_annotation, args.infer_genes, args.infer_transcripts, args.merge_strategy, args.id_spec, args.force, args.verbose, auto_convert_gtf)

    t3 = time.process_time()
    ################################
    # Step 2: Get all reference features to liftover
    ################################
    # Gene-like lift — DEFAULT since Iteration 12 (was opt-in --lift-gene-like in
    # Iteration 5). When the user did NOT supply an explicit -f/--features list,
    # auto-detect every gene-like top-level parent type (pseudogenes, ncRNA_gene,
    # structured mobile elements, ...) from the reference and write it to a temp
    # feature_types file. Setting args.features here makes BOTH this Step-2 call
    # and the vendored Liftoff invocation (Step 4, which reads args.features) lift
    # the same expanded set. --gene-only opts out (args.features stays None ->
    # get_parent_features_to_lift returns ["gene"], the pre-Iteration-12 default).
    # For a gene-only reference the auto-detect returns ["gene"], so the injected
    # file is a no-op equivalent to None (byte-identical; 24-cell matrix green).
    # --lift-gene-like is a kept no-op alias (the behaviour is now default).
    _lift_gene_like_injected = False
    if not getattr(args, "gene_only", False) and args.features is None:
        gene_like = lifton_utils.get_gene_like_feature_types(ref_db)
        auto_path = os.path.join(intermediate_dir, "auto_feature_types.txt")
        with open(auto_path, "w") as fw:
            fw.write("\n".join(gene_like) + "\n")
        args.features = auto_path
        _lift_gene_like_injected = True
        logger.log_info(f">> gene-like lift (default): auto-detected gene-like "
                        f"feature types {gene_like} (written to {auto_path}); "
                        f"pass --gene-only to lift only `gene`")
    features = lifton_utils.get_parent_features_to_lift(args.features)
    ref_features_dict, ref_features_len_dict, ref_features_reverse_dict, ref_trans_exon_num_dict = lifton_utils.get_ref_liffover_features(features, ref_db, intermediate_dir, args)

    t4 = time.process_time()
    ################################
    # Step 3: Extract protein & DNA dictionaries from the selected reference features
    ################################
    ref_trans_file = args.transcripts    
    ref_proteins_file = args.proteins    
    if (ref_proteins_file is None) or (not os.path.exists(ref_proteins_file)) or (ref_trans_file is None) or (not os.path.exists(ref_trans_file)):
        logger.log(">> Creating transcript DNA dictionary from the reference annotation ...", debug=True)
        logger.log(">> Creating transcript protein dictionary from the reference annotation ...", debug=True)
        # Phase 15b (V3.1) — streaming extractor writes FASTA directly,
        # no in-memory dict materialisation. Then re-open via pyfaidx
        # so downstream consumers see the same lazy mmap-backed
        # interface as the user-supplied -P / -T branch below.
        ref_trans_file, ref_proteins_file = \
            extract_sequence.extract_features_to_fasta(
                ref_db, features, ref_fai, intermediate_dir,
            )
        ref_trans = Fasta(ref_trans_file)
        ref_proteins = Fasta(ref_proteins_file)
    else:
        logger.log(">> Reading transcript DNA dictionary from the reference fasta ...", debug=True)
        logger.log(">> Reading transcript protein dictionary from the reference fasta ...", debug=True)
        ref_trans = Fasta(ref_trans_file)
        ref_proteins = Fasta(ref_proteins_file)
    logger.log("\t * number of transcripts: ", len(ref_trans.keys()), debug=True)
    logger.log("\t * number of proteins: ", len(ref_proteins.keys()), debug=True)
    trunc_ref_proteins = lifton_utils.get_truncated_protein(ref_proteins)
    logger.log("\t\t * number of truncated proteins: ", len(trunc_ref_proteins.keys()), debug=True)

    ################################
    # optional Step: Evaluation mode
    ################################
    if args.evaluation:
        tgt_annotation = args.output
        ref_annotation = args.reference_annotation
        print("Run LiftOn in evaluation mode")
        print("lifton_outdir     : ", lifton_outdir)
        print("Ref genome        : ", ref_genome)
        print("Target genome     : ", tgt_genome)
        print("Ref annotation    : ", args.reference_annotation)
        print("Target annotation : ", args.output)
        print("ref_trans_file    : ", ref_trans_file)
        print("ref_proteins_file : ", ref_proteins_file)
        logger.log(">> Creating target database : ", tgt_annotation, debug=True)
        os.makedirs(lifton_outdir, exist_ok=True)
        auto_convert_gtf = not args.no_auto_convert_gtf
        tgt_feature_db = annotation.Annotation(tgt_annotation, args.infer_genes, args.infer_transcripts, args.merge_strategy, args.id_spec, args.force, args.verbose, auto_convert_gtf).db_connection
        fw_score = open(lifton_outdir+"/eval.txt", "w")
        tree_dict = {}
        processed_features = 0
        for feature in features:
            for locus in tgt_feature_db.features_of_type(feature):#, limit=("chr1", 146652669, 146708545)):
                lifton_gene = run_evaluation.evaluation(None, locus, ref_db.db_connection, tgt_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, args, ENTRY_FEATURE=True)
                if processed_features % 20 == 0:
                    sys.stdout.write("\r>> LiftOn evaluated: %i features." % processed_features)
                processed_features += 1
        fw_score.close()
        return
    
    ################################
    # Step 4: Run liftoff & miniprot
    ################################
    t5 = time.process_time()
    # Output-neutral perf probe: Step 4 is subprocess-dominated, so
    # process_time (parent-CPU only) cannot measure it. Capture WALL time
    # around the whole block and emit one stderr line when LIFTON_PERF_STEP4
    # is set — lets the --parallel-aligners A/B isolate the overlap saving
    # from Step-7 noise without touching the output GFF3 or time.txt.
    _w4_start = time.perf_counter()
    if not getattr(args, "serial_aligners", False):
        # Iteration 6 (PROMOTED to default): overlap the two independent
        # external aligners so wall-clock = max(t_liftoff, t_miniprot) instead
        # of the sum. They read the same inputs and write disjoint output dirs
        # (liftoff/ vs miniprot/), consumed separately at Step 5, so the output
        # bytes are unchanged — only miniprot's *timing* moves. This is a
        # byte-neutral default flip; --serial-aligners restores the old
        # sequential path (and --parallel-aligners is a kept no-op alias).
        #
        # Liftoff MUST stay on THIS (main) thread: it reads the reference
        # gffutils DB whose SQLite connection is bound to the thread that
        # created it (Step 1), and sqlite3 forbids cross-thread use
        # (`SQLite objects created in a thread can only be used in that same
        # thread`). miniprot is an independent subprocess that never touches
        # ref_db, so IT is the one dispatched to a background worker.
        # .result() re-raises the worker's exception on the main thread,
        # preserving the serial fail-fast (run_miniprot returns None on
        # failure; run_liftoff sys.exit(1) still fires inline here).
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=1, thread_name_prefix="lifton-miniprot") as _ex:
            _f_mini = _ex.submit(
                lifton_utils.exec_miniprot, lifton_outdir, args,
                tgt_genome, ref_proteins_file)
            liftoff_annotation = lifton_utils.exec_liftoff(
                lifton_outdir, ref_db, args)
            miniprot_annotation = _f_mini.result()
        # The two tools overlapped, so the per-tool t5->t6 / t6->t7 split is
        # meaningless; collapse t6 onto t7 so the --measure_time report's
        # "Run liftoff & miniprot" line (t6 - t5) covers the whole concurrent
        # region and "Create liftoff database" (t8 - t6) stays DB-only. Under
        # this flag that entry is overlapped CPU time (process_time was never
        # a wall-clock proxy; the true wall delta lives in the A/B harness).
        t6 = time.process_time()
        t7 = t6
    else:
        liftoff_annotation = lifton_utils.exec_liftoff(lifton_outdir, ref_db, args)
        t6 = time.process_time()
        miniprot_annotation = lifton_utils.exec_miniprot(lifton_outdir, args, tgt_genome, ref_proteins_file)
        t7 = time.process_time()
    if os.environ.get("LIFTON_PERF_STEP4"):
        _mode = "serial" if getattr(args, "serial_aligners", False) else "parallel"
        sys.stderr.write(
            f"[LiftOn][perf] Step4 wall ({_mode}): "
            f"{time.perf_counter() - _w4_start:.2f}s\n")
        sys.stderr.flush()
    ################################
    # Step 5: Create liftoff and miniprot database
    ################################
    logger.log(f"\n>> Creating liftoff annotation database : {_describe_annotation_source(liftoff_annotation)}", debug=True)
    l_feature_db = annotation.Annotation(
        liftoff_annotation,
        infer_genes=False,
        infer_transcripts=False,
        merge_strategy=args.merge_strategy,
        id_spec=None,
        force=args.force,
        verbose=args.verbose,
        auto_convert_gtf=False
    ).db_connection
    t8 = time.process_time()
    logger.log(f">> Creating miniprot annotation database : {_describe_annotation_source(miniprot_annotation)}", debug=True)
    if miniprot_annotation is not None:
        m_feature_db = annotation.Annotation(
            miniprot_annotation,
            infer_genes=False,
            infer_transcripts=False,
            merge_strategy=args.merge_strategy,
            id_spec=None,
            force=args.force,
            verbose=args.verbose,
            auto_convert_gtf=False
        ).db_connection
    else:
        print(
            "[LiftOn] Skipping miniprot annotation database: miniprot produced no output.",
            file=sys.stderr,
        )
        m_feature_db = None
    fw = open(args.output, "w")
    # Phase 15a (V5.7) — emit directive prologue BEFORE any feature row.
    # Single source of truth: lifton.io.gff3_writer.format_directives.
    # Runs on the parent thread before any worker exists, so no
    # interleaving risk and no I/O lock needed.
    from lifton.io import gff3_writer as _gff3_writer
    fw.write(_gff3_writer.format_directives(
        getattr(ref_db, "directives", []) or []
    ))
    fw_score = open(f"{lifton_outdir}/score.txt", "w")
    fw_unmapped = open(f"{stats_dir}/unmapped_features.txt", "w")
    fw_extra_copy = open(f"{stats_dir}/extra_copy_features.txt", "w")
    fw_mapped_feature = open(f'{stats_dir}/mapped_feature.txt', 'w')
    fw_mapped_trans = open(f'{stats_dir}/mapped_transcript.txt', 'w')
    fw_feature_type = open(f'{stats_dir}/completeness_by_feature_type.txt', 'w')
    fw_chain = open(f"{lifton_outdir}/chain.txt", "w") if args.write_chains else None

    t9 = time.process_time()
    ################################
    # Step 6: Creating miniprot 2 Liftoff ID mapping & Initializing intervaltree
    ################################
    if m_feature_db is not None:
        ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict = lifton_utils.miniprot_id_mapping(m_feature_db)
    else:
        ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict = {}, {}
    # The interval tree gates miniprot protein-rescue suppression (Step 8): a
    # miniprot mRNA overlapping a tree locus is dropped as redundant. That set
    # must stay the protein-coding "primary" lift (`gene`) — the default gene-like
    # lift (Iter 12) adds pseudogenes/ncRNA_genes to `features` for PROCESSING,
    # but they must NOT suppress a miniprot coding-gene rescue just by overlapping
    # it (that lost CG3303/Plac8l1 etc. in the A/B). When the gene-like lift
    # injected an expanded set, seed the tree with ["gene"]; under --gene-only or
    # an explicit -f the injection didn't fire so tree == features. For a
    # gene-only ref the injected set IS ["gene"], so both paths agree
    # (byte-identical default — 24-cell matrix green).
    tree_features = ["gene"] if _lift_gene_like_injected else features
    tree_dict = intervals.initialize_interval_tree(l_feature_db, tree_features)
    transcripts_stats_dict = {'coding': {}, 'non-coding': {}, 'other': {}}
    processed_features = 0
    
    t10 = time.process_time()
    ################################
    # Step 7: Process Liftoff genes & transcripts
    #     structure 1: gene -> transcript -> exon
    #     structure 2: transcript -> exon
    #
    # Phase 9: dispatch per-locus work either serially (default) or
    # through a ThreadPoolExecutor with deterministic ordered-writer
    # buffer. The legacy serial loop is preserved byte-for-byte so the
    # default path cannot regress; the parallel path is opt-in via
    # --locus-pipeline + --threads N and emits in submission order so
    # output is byte-identical to --threads 1.
    ################################
    from lifton import parallel as _parallel
    from lifton.locus_pipeline import StepContext as _StepContext
    _ctx = _StepContext(
        ref_db=ref_db.db_connection,
        l_feature_db=l_feature_db,
        m_feature_db=m_feature_db,
        ref_id_2_m_id_trans_dict=ref_id_2_m_id_trans_dict,
        tree_dict=tree_dict,
        tgt_fai=tgt_fai,
        ref_proteins=ref_proteins,
        ref_trans=ref_trans,
        ref_features_dict=ref_features_dict,
        fw_score=fw_score,
        fw_chain=fw_chain,
        args=args,
    )
    _threads = int(getattr(args, "threads", 1) or 1)
    _use_pool = bool(getattr(args, "locus_pipeline", False)) and _threads > 1
    # Output-neutral perf probe: Step 7 is the per-locus wall-clock hot
    # spot, but it runs inside the same process so process_time mixes it
    # with worker CPU. Capture WALL time around the dispatch and emit one
    # stderr line when LIFTON_PERF_STEP7 is set — lets the Iteration-8
    # fresh-parallel A/B isolate the Step-7 speedup (no --native needed)
    # without touching the output GFF3 or time.txt.
    _w7_start = time.perf_counter()
    processed_features = _parallel.parallel_step7(
        features, l_feature_db, _ctx, fw, transcripts_stats_dict,
        threads=_threads if _use_pool else 1,
    )
    if os.environ.get("LIFTON_PERF_STEP7"):
        _mode = f"pool x{_threads}" if _use_pool else "serial"
        sys.stderr.write(
            f"[LiftOn][perf] Step7 wall ({_mode}): "
            f"{time.perf_counter() - _w7_start:.2f}s\n")
        sys.stderr.flush()

    t11 = time.process_time()
    ################################
    # Step 8: Process miniprot transcripts
    ################################
    if m_feature_db is not None:
        for mtrans in m_feature_db.features_of_type('mRNA'):
            try:
                lifton_gene = run_miniprot.process_miniprot(mtrans, ref_db, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, m_id_2_ref_id_trans_dict, ref_features_len_dict, ref_trans_exon_num_dict, ref_features_reverse_dict, args)
                if lifton_gene is None or lifton_gene.ref_gene_id is None:
                    continue
                lifton_gene.write_entry(fw, transcripts_stats_dict)
            except Exception as e:
                logger.log_error(f"Error during miniprot text output serialization ({mtrans.id}): {e}")
                
            if processed_features % 20 == 0:
                sys.stdout.write("\r>> LiftOn processed: %i features." % processed_features)
            processed_features += 1
    
    t12 = time.process_time()
    ################################
    # Step 9: Printing stats
    ################################
    try:
        stats.print_report(ref_features_dict, transcripts_stats_dict, fw_unmapped, fw_extra_copy, fw_mapped_feature, fw_mapped_trans, debug=args.debug, fw_feature_type=fw_feature_type)
    except Exception as e:
        logger.log_error(f"Failed to print report: {e}")
    finally:
        # Guarantee closure of file objects
        fw.close()
        fw_score.close()
        fw_unmapped.close()
        fw_extra_copy.close()
        fw_mapped_feature.close()
        fw_mapped_trans.close()
        fw_feature_type.close()
        if args.write_chains: fw_chain.close()

    ################################
    # Step 10: Validate output GFF3
    ################################
    if getattr(args, 'validate_output', False):
        print("\n\n*********************************************", file=sys.stderr)
        print("** Validating output GFF3                  **", file=sys.stderr)
        print("*********************************************", file=sys.stderr)
        val_result = gff3_validator.validate_gff3_file(
            gff3_path=args.output,
            check_hierarchy=True,
            check_cds_phase=True,
            check_containment=True,
            check_lifton_attrs=True,
        )
        verbose = getattr(args, 'validate_verbose', False)
        gff3_validator.print_validation_report(val_result, verbose=verbose)
        if not val_result.is_valid:
            print(
                f"\n[LiftOn] Output GFF3 has {len(val_result.errors)} error(s). "
                f"See report above for details.",
                file=sys.stderr,
            )

    t13 = time.process_time()

    if args.measure_time:
        reading_target_reference_genomes= t2 - t1
        creating_reference_annotation_database = t3 - t2
        get_all_reference_features_to_lift = t4 - t3
        extract_protein_dna_dictionaries = t5 - t4
        run_liftoff_miniprot = t6 - t5
        create_liftoff_database = t8 - t6
        create_miniprot_database = t9 - t8
        miniprot_2_liftoff_id_mapping = t10 - t9
        process_liftoff_genes_transcripts = t11 - t10
        process_miniprot_transcripts = t12 - t11
        report_stats = t13 - t12
        overall_time = t13 - t1
        print("Time taken for each step:")
        print(f"Reading target & reference genomes: {reading_target_reference_genomes}")
        print(f"Creating reference annotation database: {creating_reference_annotation_database}")
        print(f"Get all reference features to liftover: {get_all_reference_features_to_lift}")
        print(f"Extract protein & DNA dictionaries: {extract_protein_dna_dictionaries}")
        print(f"Run liftoff & miniprot: {run_liftoff_miniprot}")
        print(f"Create liftoff database: {create_liftoff_database}")
        print(f"Create miniprot database: {create_miniprot_database}")
        print(f"Miniprot 2 Liftoff ID mapping: {miniprot_2_liftoff_id_mapping}")
        print(f"Process Liftoff genes & transcripts: {process_liftoff_genes_transcripts}")
        print(f"Process miniprot transcripts: {process_miniprot_transcripts}")
        print(f"Report stats: {report_stats}")
        print(f"Overall time: {overall_time}")
        fw_time = open(f"{outdir}/time.txt", "w")
        fw_time.write(f"{reading_target_reference_genomes}\tReading target & reference genomes\n")
        fw_time.write(f"{creating_reference_annotation_database}\tCreating reference annotation database\n")
        fw_time.write(f"{get_all_reference_features_to_lift}\tGet all reference features to liftover\n")
        fw_time.write(f"{extract_protein_dna_dictionaries}\tExtract protein & DNA dictionaries\n")
        fw_time.write(f"{run_liftoff_miniprot}\tRun liftoff & miniprot\n")
        fw_time.write(f"{create_liftoff_database}\tCreate liftoff database\n")
        fw_time.write(f"{create_miniprot_database}\tCreate miniprot database\n")
        fw_time.write(f"{miniprot_2_liftoff_id_mapping}\tMiniprot 2 Liftoff ID mapping\n")
        fw_time.write(f"{process_liftoff_genes_transcripts}\tProcess Liftoff genes & transcripts\n")
        fw_time.write(f"{process_miniprot_transcripts}\tProcess miniprot transcripts\n")
        fw_time.write(f"{report_stats}\tReport stats\n")
        fw_time.write(f"{overall_time}\tOverall time\n")
        fw_time.close()


def main(arglist=None):
    banner = '''
====================================================================
An accurate homology lift-over tool between assemblies
====================================================================


    ██╗     ██╗███████╗████████╗ ██████╗ ███╗   ██╗
    ██║     ██║██╔════╝╚══██╔══╝██╔═══██╗████╗  ██║
    ██║     ██║█████╗     ██║   ██║   ██║██╔██╗ ██║
    ██║     ██║██╔══╝     ██║   ██║   ██║██║╚██╗██║
    ███████╗██║██║        ██║   ╚██████╔╝██║ ╚████║
    ╚══════╝╚═╝╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═══╝
    '''
    print(banner, file=sys.stderr)
    args = parse_args(arglist)
    if not run_miniprot.check_miniprot_installed(): 
        sys.exit("miniprot is not installed. Please install miniprot before running LiftOn.")
    run_all_lifton_steps(args)
