from lifton import intervals, lifton_utils, annotation, extract_sequence, stats, logger, run_liftoff, run_miniprot, run_evaluation, gff3_validator, __version__
from intervaltree import IntervalTree
import argparse
from pyfaidx import Fasta
import os, sys
import re
import time
import concurrent.futures
from contextlib import redirect_stdout
from urllib.parse import unquote

from lifton.exceptions import (
    LiftOnError,
    LiftOnPartialOutputError,
    LiftOnValidationError,
)
from lifton.run_manifest import RunManifest
from lifton.locus_pipeline import safe_exception_text


def _allocator_source_ids(database):
    """Collect only declared IDs relevant to file-scoped CDS allocation.

    Direct database and in-memory annotation inputs do not have an
    ``AnnotationScanResult``. Inspect logical ``ID`` attributes rather than
    database row keys: merge strategies may suffix internal keys for valid
    discontinuous CDS rows.
    """

    cds_namespace_ids = set()
    copy_suffix_ids = set()
    try:
        features = database.all_features()
    except Exception:
        return cds_namespace_ids, copy_suffix_ids
    for feature in features:
        values = getattr(feature, "attributes", {}).get("ID", [])
        if isinstance(values, str):
            values = [values]
        for value in values:
            identifier = unquote(str(value))
            if identifier.startswith("cds-"):
                cds_namespace_ids.add(identifier)
            if re.search(r"_\d+$", identifier):
                copy_suffix_ids.add(identifier)
    return cds_namespace_ids, copy_suffix_ids


def _run_outdirs(output):
    """Return the user-output directory and LiftOn artifact directory."""
    if output == "stdout":
        outdir = "."
    else:
        outdir = os.path.dirname(output) or "."
    return outdir, os.path.join(outdir, "lifton_output")


def _console_stream(args):
    """Keep stdout machine-readable when it is the GFF3 destination."""
    return sys.stderr if getattr(args, "output", None) == "stdout" else sys.stdout


def _ensure_run_manifest(args, outdir, lifton_outdir):
    """Create one provenance collector per run and attach it to ``args``."""
    manifest = getattr(args, "_run_manifest", None)
    if manifest is not None:
        return manifest

    inputs = {
        "target_genome": getattr(args, "target", None),
        "reference_genome": getattr(args, "reference", None),
        "reference_annotation": getattr(args, "reference_annotation", None),
        "reference_proteins": getattr(args, "proteins", None),
        "reference_transcripts": getattr(args, "transcripts", None),
        "precomputed_liftoff": getattr(args, "liftoff", None),
        "precomputed_miniprot": getattr(args, "miniprot", None),
    }
    minimap2 = getattr(args, "m", None) or "minimap2"
    manifest = RunManifest(
        argv=getattr(args, "_argv", sys.argv),
        options=vars(args),
        inputs=inputs,
        backend={
            "stream": bool(getattr(args, "stream", False)),
            "inmemory_liftoff": bool(
                getattr(args, "inmemory_liftoff", False)
            ),
            "native": bool(getattr(args, "native", False)),
            "locus_pipeline": bool(getattr(args, "locus_pipeline", False)),
        },
        tool_commands={
            "minimap2": (minimap2, "--version"),
            "miniprot": ("miniprot", "--version"),
        },
        git_cwd=os.getcwd(),
    )
    args._run_manifest = manifest
    args._run_manifest_path = os.path.join(
        lifton_outdir, "run_manifest.json"
    )
    args._active_manifest_phase = None
    args._pipeline_failures = []
    return manifest


def _switch_manifest_phase(args, name, details=None):
    manifest = getattr(args, "_run_manifest", None)
    if manifest is None:
        return
    active = getattr(args, "_active_manifest_phase", None)
    if active is not None:
        manifest.end_phase(active)
    args._active_manifest_phase = name
    if name is not None:
        manifest.start_phase(name, details)


def _record_pipeline_failure(args, phase, error, *, details=None,
                             fatal=False, affects_completeness=True):
    details = dict(details or {})
    message = safe_exception_text(error, details.get("traceback"))
    record = {
        "phase": str(phase),
        "message": message,
        "details": details,
        "fatal": bool(fatal),
    }
    if affects_completeness:
        getattr(args, "_pipeline_failures", []).append(record)
    manifest = getattr(args, "_run_manifest", None)
    if manifest is not None:
        manifest.record_failure(
            phase, message, fatal=fatal, details=record["details"]
        )


def _finish_run_manifest(args, status):
    manifest = getattr(args, "_run_manifest", None)
    path = getattr(args, "_run_manifest_path", None)
    if manifest is None or path is None:
        return
    active = getattr(args, "_active_manifest_phase", None)
    if active is not None:
        try:
            manifest.end_phase(active, status=(
                "success" if status in {"success", "partial_success"}
                else "failed"
            ))
        except ValueError:
            pass
        args._active_manifest_phase = None
    manifest.finish(status)
    manifest.write(path)


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
        help='The maximum length ratio of a protein-coding transcript to the longest protein-coding transcript within a gene locus, as identified exclusively by miniprot in the target genome, is set by default to MAX_MINIPROT=1.5.'
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
        '--allow-partial-output', action='store_true', default=False,
        help='Publish an output even when one or more loci could not be '
             'serialized. The default reports partial results as a failed run; '
             'all skipped loci are recorded in run_manifest.json.'
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
        help='Locus-major fan-out: dispatch per-locus Step 7, Step 8, and '
             'evaluation work through bounded workers sized by --threads. '
             'Output is emitted in submission order so --threads N is '
             'byte-identical to --threads 1; this flag changes scheduling, '
             'not algorithms. As of Iteration 8 this works on the DEFAULT '
             '(gffutils) backend WITHOUT --native — per-locus work runs '
             'against materialised proxy DBs, so any backend is thread-safe '
             '(in-memory gffbase databases use independent DuckDB cursors; '
             'other non-reopenable databases are materialised on the parent '
             'thread before worker dispatch). '
             'Set LIFTON_PARALLEL_BLOCK_GFFUTILS=1 to opt back out to '
             'serial-on-gffutils. Note: combined with the default '
             'concurrent Step 4, peak busy cores can reach ~N+1.'
    )
    parser.add_argument(
        '--step7-max-inflight', dest='step7_max_inflight', type=int,
        default=None, metavar='N',
        help='Bound submitted-but-not-emitted Step-7 loci. The default is '
             '2 * --threads; lower values reduce peak memory and higher values '
             'may improve utilization for uneven loci.'
    )
    parser.add_argument(
        '--step8-max-inflight', dest='step8_max_inflight', type=int,
        default=None, metavar='N',
        help='Bound analyzed-but-not-published miniprot candidates in Step 8. '
             'The default is 2 * --threads; lower values reduce peak memory.'
    )
    parser.add_argument(
        '--evaluation-max-inflight', dest='evaluation_max_inflight', type=int,
        default=None, metavar='N',
        help='Bound submitted-but-not-written loci in evaluation mode. The '
             'default is 2 * --threads; lower values reduce peak memory.'
    )
    parser.add_argument(
        '--native', dest='native', action='store_true', default=False,
        help='Enable experimental native compatibility hooks. The proven '
             'miniprot subprocess/direct-stream and minimap2 paths remain '
             'the defaults, and bounded locus workers do not require this '
             'flag. Set LIFTON_NATIVE_LIFTOFF_ALIGN=1 as well to opt into '
             'the mappy Liftoff path (for example when minimap2 is absent); '
             'it falls back gracefully when mappy is unavailable.'
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
        '--no-miniprot-rescue', dest='miniprot_rescue', action='store_false', default=True,
        help='Opt OUT of the miniprot-only rescue pass. As of the Iteration-23 '
             'promotion the rescue is the DEFAULT: for a reference coding gene '
             'the DNA lift MISSED ENTIRELY (its miniprot mRNA overlaps no lifted '
             'gene locus), LiftOn emits the miniprot model even when its length '
             'ratio falls outside the default -min_miniprot/-max_miniprot band -- '
             'gated instead by a protein-identity floor '
             '(LIFTON_MINIPROT_RESCUE_MIN_ID, default 0.5) within a wider sanity '
             'band (LIFTON_MINIPROT_RESCUE_LEN=lo,hi, default 0.5,2.0). The '
             'rescue runs as a SEPARATE pass AFTER Step 8 closes, with a '
             'ref-gene-id dedup set, so the default Step-7+8 output is '
             'byte-identical between OFF and ON (off is a strict subset of on) '
             'and 0-redundant; it recovers genuinely-missing genes at large '
             'evolutionary distance (net +recall on the distant/very-distant '
             'tier -- see benchmarks/compare/miniprot_rescue_ab.md, 8/8 datasets '
             '0 lost / 0 redundant). Pass --no-miniprot-rescue (or env '
             'LIFTON_MINIPROT_RESCUE=0) to restore the pre-Iteration-23 lift '
             'with no miniprot-only rescue.'
    )
    parser.add_argument(
        '--miniprot-rescue', dest='miniprot_rescue_alias', action='store_true', default=False,
        help='No-op alias (kept for backward compatibility). The miniprot-only '
             'rescue is now the DEFAULT (Iteration-23 promotion), so '
             '--miniprot-rescue has no effect; pass --no-miniprot-rescue to opt '
             'OUT. Env LIFTON_MINIPROT_RESCUE=1 force-enables, =0 force-disables.'
    )
    parser.add_argument(
        '--miniprot-cross-locus-rescue', dest='cross_locus_rescue',
        action='store_true', default=False,
        help='[EXPERIMENTAL, opt-in] Phase D cross-locus rescue. For a reference '
             'coding gene whose DNA lift produced only a WEAK model (best emitted '
             'protein identity < LIFTON_CROSS_LOCUS_MAX_LIFTOFF, default 0.5) '
             'while miniprot aligns it cleanly to a DIFFERENT chromosome '
             '(paralog/duplicate disagreement -- the Figure-4 very-distant '
             'residual candidate-3 cannot reach), REPLACE the weak gene: drop all '
             'its blocks from the output GFF3 and append the miniprot model '
             '(tagged lifton_rescue=cross_locus) iff its identity beats the '
             'emitted best by >= LIFTON_CROSS_LOCUS_MIN_GAIN (default 0.3). '
             'Replace-not-add => duplicate-safe. Runs as a separate post-write '
             'pass, so OFF (the default) is byte-identical. Env '
             'LIFTON_CROSS_LOCUS_RESCUE=1 force-enables. SYNTENY CAVEAT: trades '
             'synteny for protein identity -- keep opt-in until a strict A/B proof.'
    )
    parser.add_argument(
        '--no-miniprot-candidate', dest='miniprot_candidate',
        action='store_false', default=True,
        help='Opt OUT of the miniprot-only best-of-outcome candidate (the 3rd '
             'candidate, DEFAULT-ON). For a transcript where the DNA lift is '
             'imperfect and miniprot overlaps the SAME locus, LiftOn also scores '
             "miniprot's NATIVE model (a clean CDS-only scaffold) and keeps it "
             'ONLY when its ORF-rescued protein identity is STRICTLY better than '
             'the 2-way winner max(Liftoff+ORF, chained-merge+ORF) -- so '
             'per-transcript identity is non-decreasing (additive) and close '
             'pairs barely change. It recovers the very-distant residual where '
             'the DNA lift collapses to a truncated stub while miniprot aligns '
             'the protein near-perfectly (see benchmarks/compare/'
             'miniprot_candidate_ab.md). Pass --no-miniprot-candidate (or env '
             'LIFTON_MINIPROT_CANDIDATE=0) to restore the 2-way merge.'
    )
    parser.add_argument(
        '--miniprot-candidate', dest='miniprot_candidate_alias',
        action='store_true', default=False,
        help='No-op alias (kept for backward compatibility). The miniprot-only '
             'candidate is now the DEFAULT, so --miniprot-candidate has no '
             'effect; pass --no-miniprot-candidate to opt OUT. Env '
             'LIFTON_MINIPROT_CANDIDATE=1 force-enables, =0 force-disables.'
    )
    parser.add_argument(
        '--no-adaptive-rescue-floor', dest='adaptive_rescue_floor',
        action='store_false', default=True,
        help='Opt OUT of the divergence-adaptive miniprot-only-rescue floor '
             '(DEFAULT-ON). The rescue admits a miniprot-only gene when its '
             'protein identity clears a floor; by default that floor is LOWERED '
             'toward 0.30 as the DNA-lift gene recall drops (divergent pairs) to '
             'recover more genuinely-missing genes, and stays at the base '
             '(LIFTON_MINIPROT_RESCUE_MIN_ID, default 0.5) on same/close-species '
             '(high recall). Pass --no-adaptive-rescue-floor (or env '
             'LIFTON_RESCUE_ADAPTIVE_FLOOR=0) to restore the FIXED base floor. '
             'Thresholds: LIFTON_RESCUE_FLOOR_MIN/R_LOW/R_HIGH.'
    )
    parser.add_argument(
        '--adaptive-rescue-floor', dest='adaptive_rescue_floor_alias',
        action='store_true', default=False,
        help='No-op alias (kept for backward compatibility). The adaptive rescue '
             'floor is now the DEFAULT, so --adaptive-rescue-floor has no effect; '
             'pass --no-adaptive-rescue-floor to opt OUT. Env '
             'LIFTON_RESCUE_ADAPTIVE_FLOOR=1 force-enables, =0 force-disables.'
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
    args._argv = ([sys.argv[0], *sys.argv[1:]] if arglist is None
                  else ["lifton", *list(arglist)])
    if args.evaluation_liftoff_chm13:
        args.evaluation = True
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
    if args.step7_max_inflight is not None and args.step7_max_inflight < 1:
        parser.error("--step7-max-inflight must be at least 1")
    if args.step8_max_inflight is not None and args.step8_max_inflight < 1:
        parser.error("--step8-max-inflight must be at least 1")
    if (args.evaluation_max_inflight is not None
            and args.evaluation_max_inflight < 1):
        parser.error("--evaluation-max-inflight must be at least 1")
    return args
    

def resolve_miniprot_rescue_args(args):
    """Resolve the miniprot-only-rescue flag + its tunables ONCE, honouring the
    env escape hatches, and stash them on `args` so the post-Step-8 rescue pass
    (lifton.miniprot_rescue) can read them. Pure/idempotent (no I/O beyond env
    reads) so it is unit-testable.

    As of the Iteration-23 promotion the rescue is DEFAULT-ON (a separate
    post-Step-8 pass => default Step-7+8 output is byte-identical OFF-vs-ON).
    Opt out with --no-miniprot-rescue (=> args.miniprot_rescue False) or
    LIFTON_MINIPROT_RESCUE=0.

    - LIFTON_MINIPROT_RESCUE=1|true|yes   -> force ON
    - LIFTON_MINIPROT_RESCUE=0|false|no   -> force OFF
    - LIFTON_MINIPROT_RESCUE_MIN_ID=<f>   -> protein-identity floor (default 0.5)
    - LIFTON_MINIPROT_RESCUE_LEN=lo,hi    -> wider sanity band (default 0.5,2.0)
    """
    rescue = getattr(args, "miniprot_rescue", True)
    _env = os.environ.get("LIFTON_MINIPROT_RESCUE", "").lower()
    if _env in ("1", "true", "yes"):
        rescue = True
    elif _env in ("0", "false", "no"):
        rescue = False
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
    # Divergence-adaptive rescue floor (PROMOTED default ON). Env
    # LIFTON_RESCUE_ADAPTIVE_FLOOR wins, else args.adaptive_rescue_floor (default
    # True); --no-adaptive-rescue-floor sets it False (restores the fixed 0.5
    # floor). Read at apply time by miniprot_rescue._adaptive_floor_on (env wins
    # there too), so this only sets the default used when the env var is unset.
    _env_af = os.environ.get("LIFTON_RESCUE_ADAPTIVE_FLOOR")
    if _env_af is not None:
        args.adaptive_rescue_floor = _env_af.lower() not in ("0", "false", "no", "")
    else:
        args.adaptive_rescue_floor = bool(getattr(args, "adaptive_rescue_floor", True))
    return args


def resolve_miniprot_candidate_args(args):
    """Resolve the miniprot-only best-of-outcome candidate (candidate-3) flag
    ONCE and publish it to run_liftoff so the per-transcript decision point
    (process_liftoff_with_protein, which does not receive `args`) can read it via
    the module default. Pure/idempotent; unit-testable.

    Candidate-3 is DEFAULT-ON. --no-miniprot-candidate => args.miniprot_candidate
    False. The LIFTON_MINIPROT_CANDIDATE env var is honoured at READ time inside
    run_liftoff._miniprot_candidate_enabled (env WINS), so the 24-cell matrix and
    the A/B harness (which force via env) are unaffected by the CLI default; this
    resolver only sets the default used when the env var is unset.
    """
    from lifton import run_liftoff as _run_liftoff
    enabled = getattr(args, "miniprot_candidate", True)
    args.miniprot_candidate = enabled
    _run_liftoff._MINIPROT_CANDIDATE_DEFAULT = enabled
    return args


def run_all_lifton_steps(args):
    t1 = time.process_time()
    # Iteration-3 "band everything" alignment is the DEFAULT (set at align-module
    # import). --full-dp-align restores the exact pre-Iteration-3 giant-only path
    # for manuscript reproduction / paranoid accuracy. Lazy import keeps `align`
    # (which pulls in the lifton_class↔lifton_utils cycle) off the module top.
    # Configure both arms explicitly so repeated in-process calls cannot leak a
    # previous run's ``--full-dp-align`` setting into the next invocation.
    from lifton import align as _align
    _align.configure_alignment(band=not getattr(args, "full_dp_align", False))
    resolve_miniprot_rescue_args(args)
    resolve_miniprot_candidate_args(args)
    ################################
    # Step 0: Reading target & reference genomes
    ################################
    tgt_genome = args.target
    ref_genome = args.reference
    outdir, lifton_outdir = _run_outdirs(args.output)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(lifton_outdir, exist_ok=True)
    args.directory = "intermediate_files/"
    intermediate_dir = f"{outdir}/lifton_output/{args.directory}"
    os.makedirs(intermediate_dir, exist_ok=True)
    stats_dir= f"{outdir}/lifton_output/stats/"
    os.makedirs(stats_dir, exist_ok=True)
    args.directory = intermediate_dir
    manifest = _ensure_run_manifest(args, outdir, lifton_outdir)
    _switch_manifest_phase(args, "read_genomes")
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
    _switch_manifest_phase(args, "validate_reference_annotation")
    reference_seqids = set(ref_fai.keys())
    reference_annotation_scan = None
    if annotation.is_annotation_database_file(args.reference_annotation):
        findings = []
        manifest.set_backend_choice("reference_annotation_input", "database")
    else:
        reference_annotation_scan = annotation.scan_annotation(
            args.reference_annotation,
            target_seqids=reference_seqids,
        )
        findings = list(reference_annotation_scan.ncbi_findings)
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
        _record_pipeline_failure(
            args, "validate_reference_annotation",
            "reference annotation failed strict GFF3 validation",
            details={"errors": sum(
                1 for finding in findings if finding.severity == "error"
            )},
            fatal=True,
        )
        _finish_run_manifest(args, "failed")
        sys.exit(2)
    ################################
    # Step 1: Building database from the reference annotation
    ################################
    _switch_manifest_phase(args, "build_reference_database")
    logger.log("\n>> Creating reference annotation database : ", args.reference_annotation, debug=True)
    auto_convert_gtf = not args.no_auto_convert_gtf
    ref_db = annotation.Annotation(
        args.reference_annotation,
        args.infer_genes,
        args.infer_transcripts,
        args.merge_strategy,
        args.id_spec,
        args.force,
        args.verbose,
        auto_convert_gtf,
        scan_result=reference_annotation_scan,
    )
    manifest.set_backend_choice("reference_annotation", ref_db.backend)
    manifest.set_cache_choice(
        "reference_annotation", getattr(ref_db, "cache_status", "unknown")
    )

    t3 = time.process_time()
    ################################
    # Step 2: Get all reference features to liftover
    ################################
    _switch_manifest_phase(args, "select_reference_features")
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
    _switch_manifest_phase(args, "load_reference_sequences")
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
    manifest.record_count("reference_transcripts", len(ref_trans.keys()))
    manifest.record_count("reference_proteins", len(ref_proteins.keys()))
    manifest.record_count("truncated_reference_proteins", len(trunc_ref_proteins.keys()))

    ################################
    # optional Step: Evaluation mode
    ################################
    if args.evaluation:
        _switch_manifest_phase(args, "evaluation")
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
        from lifton.output_transaction import OutputTransaction
        evaluation_transaction = OutputTransaction(
            os.path.join(lifton_outdir, "eval.txt")
        )
        evaluation_transaction.install_signal_handlers()
        args._output_transaction = evaluation_transaction
        fw_score = evaluation_transaction.open()
        tree_dict = {}
        processed_features = 0
        evaluation_failures = []
        evaluation_threads = int(getattr(args, "threads", 1) or 1)
        if (bool(getattr(args, "locus_pipeline", False))
                and evaluation_threads > 1):
            def evaluation_loci():
                for feature in features:
                    yield from tgt_feature_db.features_of_type(feature)

            processed_features, evaluation_failures = \
                run_evaluation.parallel_evaluate_loci(
                    evaluation_loci(),
                    ref_db.db_connection,
                    tgt_feature_db,
                    tree_dict,
                    tgt_fai,
                    ref_proteins,
                    ref_trans,
                    ref_features_dict,
                    fw_score,
                    args,
                    threads=evaluation_threads,
                    max_inflight=getattr(
                        args, "evaluation_max_inflight", None,
                    ),
                )
        else:
            for feature in features:
                for locus in tgt_feature_db.features_of_type(feature):
                    run_evaluation.evaluation(
                        None, locus, ref_db.db_connection, tgt_feature_db,
                        tree_dict, tgt_fai, ref_proteins, ref_trans,
                        ref_features_dict, fw_score, args,
                        ENTRY_FEATURE=True,
                    )
                    if processed_features % 20 == 0:
                        sys.stdout.write(
                            "\r>> LiftOn evaluated: %i features."
                            % processed_features
                        )
                    processed_features += 1
        evaluation_transaction.close()
        for failure in evaluation_failures:
            _record_pipeline_failure(
                args,
                "evaluation",
                failure.error,
                details={
                    "locus": failure.locus_id,
                    "traceback": failure.error_tb,
                },
                fatal=True,
            )
        manifest.record_count("evaluated_features", processed_features)
        manifest.record_count(
            "evaluation_failures", len(evaluation_failures),
        )
        manifest.record_count(
            "evaluation_max_inflight_observed",
            getattr(
                args,
                "_evaluation_max_inflight_observed",
                1 if processed_features else 0,
            ),
        )
        if evaluation_failures:
            partial_path = evaluation_transaction.abort()
            args._output_transaction = None
            logger.log_error(
                f"Incomplete evaluation was preserved at {partial_path}"
            )
            _switch_manifest_phase(args, None)
            _finish_run_manifest(args, "failed")
            raise LiftOnError(
                f"evaluation failed for {len(evaluation_failures)} loci"
            )
        evaluation_transaction.commit()
        args._output_transaction = None
        _switch_manifest_phase(args, None)
        _finish_run_manifest(args, "success")
        return
    
    ################################
    # Step 4: Run liftoff & miniprot
    ################################
    _switch_manifest_phase(args, "run_aligners", {
        "parallel": not bool(getattr(args, "serial_aligners", False)),
    })
    t5 = time.process_time()
    # Output-neutral perf probe: Step 4 is subprocess-dominated, so
    # process_time (parent-CPU only) cannot measure it. Capture WALL time
    # around the whole block and emit one stderr line when LIFTON_PERF_STEP4
    # is set — lets the --parallel-aligners A/B isolate the overlap saving
    # from Step-7 noise without touching the output GFF3 or time.txt.
    _w4_start = time.perf_counter()
    if len(ref_proteins.keys()) == 0:
        logger.log_info(
            ">> Reference protein set is empty; skipping miniprot."
        )
        liftoff_annotation = lifton_utils.exec_liftoff(
            lifton_outdir, ref_db, args
        )
        miniprot_annotation = None
        t6 = time.process_time()
        t7 = t6
    elif not getattr(args, "serial_aligners", False):
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
    if miniprot_annotation is None and len(ref_proteins.keys()) > 0:
        _record_pipeline_failure(
            args, "run_aligners",
            "miniprot produced no usable annotation for a non-empty protein set",
            details={"reference_proteins": len(ref_proteins.keys())},
        )
    if os.environ.get("LIFTON_PERF_STEP4"):
        _mode = "serial" if getattr(args, "serial_aligners", False) else "parallel"
        sys.stderr.write(
            f"[LiftOn][perf] Step4 wall ({_mode}): "
            f"{time.perf_counter() - _w4_start:.2f}s\n")
        sys.stderr.flush()
    ################################
    # Step 5: Create liftoff and miniprot database
    ################################
    _switch_manifest_phase(args, "build_alignment_databases")
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
    )
    manifest.set_backend_choice(
        "liftoff_annotation", l_feature_db.backend
    )
    manifest.set_cache_choice(
        "liftoff_annotation",
        getattr(l_feature_db, "cache_status", "unknown"),
    )
    l_feature_db = l_feature_db.db_connection
    t8 = time.process_time()
    logger.log(f">> Creating miniprot annotation database : {_describe_annotation_source(miniprot_annotation)}", debug=True)
    if miniprot_annotation is not None:
        m_annotation = annotation.Annotation(
            miniprot_annotation,
            infer_genes=False,
            infer_transcripts=False,
            merge_strategy=args.merge_strategy,
            id_spec=None,
            force=args.force,
            verbose=args.verbose,
            auto_convert_gtf=False
        )
        manifest.set_backend_choice("miniprot_annotation", m_annotation.backend)
        manifest.set_cache_choice(
            "miniprot_annotation",
            getattr(m_annotation, "cache_status", "unknown"),
        )
        m_feature_db = m_annotation.db_connection
    else:
        print(
            "[LiftOn] Skipping miniprot annotation database: miniprot produced no output.",
            file=sys.stderr,
        )
        m_annotation = None
        m_feature_db = None
    # Reserve every explicit source identifier in LiftOn's synthesized CDS
    # namespace before the first output row. This makes ID allocation
    # file-scoped: an early ID-less transcript cannot take a stable ID that a
    # later transcript (or copy) needs to preserve.
    from lifton.cds_id_allocator import CdsIdAllocator
    source_scan = getattr(ref_db, "scan_result", None)
    if source_scan is not None:
        reserved_cds_ids = set(source_scan.cds_namespace_ids)
        reserved_copy_suffix_ids = set(source_scan.copy_suffix_ids)
    else:
        reserved_cds_ids, reserved_copy_suffix_ids = (
            _allocator_source_ids(ref_db.db_connection)
        )

    args._cds_id_allocator = CdsIdAllocator(
        reserved_cds_ids,
        exact_reserved_source_ids=reserved_copy_suffix_ids,
    )
    manifest.record_count(
        "reserved_source_cds_namespace_ids", len(reserved_cds_ids)
    )
    manifest.record_count(
        "reserved_source_copy_suffix_ids",
        len(reserved_copy_suffix_ids),
    )
    from lifton.output_transaction import OutputTransaction
    output_transaction = OutputTransaction(
        args.output,
        stdout=getattr(args, "_gff_stdout", None),
        work_dir=lifton_outdir,
    )
    output_transaction.install_signal_handlers()
    args._output_transaction = output_transaction
    staged_output_path = os.fspath(output_transaction.working_path)
    fw = output_transaction.open()
    # Phase 15a (V5.7) — emit directive prologue BEFORE any feature row.
    # Single source of truth: lifton.io.gff3_writer.format_directives.
    # Runs on the parent thread before any worker exists, so no
    # interleaving risk and no I/O lock needed.
    from lifton.io import gff3_writer as _gff3_writer
    fw.write(_gff3_writer.format_directives(
        _gff3_writer.target_output_directives(
            getattr(ref_db, "directives", []) or []
        )
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
    _switch_manifest_phase(args, "initialize_locus_state")
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
    _switch_manifest_phase(args, "process_liftoff_loci")
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
    manifest.record_count(
        "step7_max_inflight_observed",
        getattr(args, "_step7_max_inflight_observed", 0),
    )
    step7_strategy = getattr(args, "_step7_parallel_strategy", None)
    if step7_strategy:
        manifest.set_backend_choice("step7_dispatch", step7_strategy)
    step7_materialisation_constraint = getattr(
        args, "_step7_materialisation_constraint", None,
    )
    if step7_materialisation_constraint:
        manifest.set_backend_choice(
            "step7_materialisation_constraint",
            step7_materialisation_constraint,
        )
    step7_fallback_reason = getattr(
        args, "_step7_parallel_fallback_reason", None,
    )
    if step7_fallback_reason:
        manifest.set_backend_choice(
            "step7_parallel_fallback", step7_fallback_reason,
        )
    not_emittable_count = 0
    for failure in getattr(_ctx, "failure_records", []):
        kind = failure.get("kind", "processing_error")
        if kind == "not_emittable":
            not_emittable_count += 1
            continue
        _record_pipeline_failure(
            args,
            "process_liftoff_loci",
            failure.get("message", kind),
            details={
                key: value for key, value in failure.items()
                if key not in {"message"}
            },
        )
    manifest.record_count("liftoff_loci_not_emittable", not_emittable_count)
    if os.environ.get("LIFTON_PERF_STEP7"):
        _mode = f"pool x{_threads}" if _use_pool else "serial"
        sys.stderr.write(
            f"[LiftOn][perf] Step7 wall ({_mode}): "
            f"{time.perf_counter() - _w7_start:.2f}s\n")
        sys.stderr.flush()

    # Iteration 23: clean separate-pass miniprot-only rescue (flag-gated).
    # Harvest the reference gene ids the DNA lift (Step 7) emitted, using the
    # SAME resolver Step 7 used (run_liftoff.py:111), so the post-Step-8 rescue
    # pass can skip already-lifted genes -- the dedup that makes it 0-redundant
    # + off ⊆ on. When --miniprot-rescue is OFF this set stays empty and is
    # never consulted, so the default path is byte-identical.
    emitted_ref_gene_ids = set(getattr(
        _ctx, "emitted_ref_gene_ids", set()
    ))
    # Step 7 needs a tree of every Liftoff locus when pairing DNA and protein
    # evidence. Step 8 has a different contract: suppress miniprot only where a
    # hierarchy was actually emitted. Rebuild that suppression tree from the
    # ordered commit ledger so an un-emittable Liftoff attempt cannot hide a
    # valid miniprot rescue.
    tree_dict = {}
    from lifton.intervals import _make_interval
    for seqid, start, end, gene_id in getattr(
            _ctx, "emitted_intervals", []):
        tree_dict.setdefault(seqid, IntervalTree()).add(
            _make_interval(start, end, gene_id)
        )
    manifest.record_count(
        "liftoff_genes_emitted", len(_ctx.emitted_intervals)
    )

    t11 = time.process_time()
    ################################
    # Step 8: Process miniprot transcripts
    ################################
    _switch_manifest_phase(args, "process_miniprot_loci")
    step8_candidates_processed = 0
    step8_genes_emitted = 0
    step8_max_inflight_observed = 0
    step8_pre_overlap_elided = 0
    step8_analysis_submitted = 0
    step8_child_batch_calls = 0
    step8_child_batch_fallbacks = 0
    step8_child_scalar_materializations = 0
    if m_feature_db is None:
        step8_dispatch = "skipped"
        step8_overlap_filter = "skipped"
    elif _use_pool:
        step8_dispatch = (
            "parallel_batched_parent_materialise"
            if callable(getattr(
                m_feature_db, "children_batched_features", None,
            ))
            else "parallel_scalar_parent_materialise"
        )
        step8_overlap_filter = "monotonic_precheck_ordered_recheck"
    else:
        step8_dispatch = "serial"
        step8_overlap_filter = "in_process"
    manifest.set_backend_choice("step8_dispatch", step8_dispatch)
    manifest.set_backend_choice(
        "step8_overlap_filter", step8_overlap_filter,
    )
    if m_feature_db is not None and _use_pool:
        from lifton import miniprot_pipeline as _miniprot_pipeline

        step8_outcome = _miniprot_pipeline.parallel_step8(
            m_feature_db.features_of_type("mRNA"),
            ref_db,
            m_feature_db,
            tree_dict,
            tgt_fai,
            ref_proteins,
            ref_trans,
            ref_features_dict,
            fw,
            fw_score,
            transcripts_stats_dict,
            m_id_2_ref_id_trans_dict,
            ref_features_len_dict,
            ref_trans_exon_num_dict,
            ref_features_reverse_dict,
            args,
            threads=_threads,
            max_inflight=getattr(args, "step8_max_inflight", None),
        )
        processed_features += step8_outcome.processed
        emitted_ref_gene_ids.update(step8_outcome.emitted_ref_gene_ids)
        step8_candidates_processed = step8_outcome.processed
        step8_genes_emitted = step8_outcome.emitted
        step8_max_inflight_observed = step8_outcome.max_inflight_observed
        step8_pre_overlap_elided = step8_outcome.preoverlap_elided
        step8_analysis_submitted = step8_outcome.analysis_submitted
        step8_child_batch_calls = step8_outcome.child_batch_calls
        step8_child_batch_fallbacks = step8_outcome.child_batch_fallbacks
        step8_child_scalar_materializations = (
            step8_outcome.child_scalar_materializations
        )
        for failure in step8_outcome.failures:
            _record_pipeline_failure(
                args,
                "process_miniprot_loci",
                failure.get("message", failure.get("kind", "failure")),
                details={
                    key: value for key, value in failure.items()
                    if key != "message"
                },
            )
    elif m_feature_db is not None:
        from lifton.locus_pipeline import (
            DeferredStateJournal, commit_locus_delta,
        )
        for mtrans in m_feature_db.features_of_type('mRNA'):
            step8_candidates_processed += 1
            processed_features += 1
            if processed_features % 20 == 0:
                _console_stream(args).write(
                    "\r>> LiftOn processed: %i features."
                    % processed_features
                )
            state_journal = DeferredStateJournal(
                ref_features_dict, buffer_score=True,
            )
            try:
                lifton_gene = run_miniprot.process_miniprot(
                    mtrans, ref_db, m_feature_db, tree_dict, tgt_fai,
                    ref_proteins, ref_trans, ref_features_dict,
                    state_journal.score_handle,
                    m_id_2_ref_id_trans_dict, ref_features_len_dict,
                    ref_trans_exon_num_dict, ref_features_reverse_dict, args,
                    state_journal=state_journal,
                )
                if lifton_gene is None or lifton_gene.ref_gene_id is None:
                    continue
                cds_id_allocator = getattr(
                    args, "_cds_id_allocator", None
                )
                write_result = (
                    lifton_gene.write_entry(fw, transcripts_stats_dict)
                    if cds_id_allocator is None
                    else lifton_gene.write_entry(
                        fw,
                        transcripts_stats_dict,
                        cds_id_allocator=cds_id_allocator,
                    )
                )
                serialization_failures = getattr(
                    lifton_gene, "_serialization_failures", []
                )
                for feature_id, detail in serialization_failures:
                    _record_pipeline_failure(
                        args, "process_miniprot_loci", detail,
                        details={
                            "mRNA": mtrans.id,
                            "feature_id": feature_id,
                            "kind": "serialization_error",
                        },
                    )
                if write_result is False:
                    if not serialization_failures:
                        _record_pipeline_failure(
                            args, "process_miniprot_loci",
                            "miniprot hierarchy was not serializable",
                            details={"mRNA": mtrans.id},
                        )
                    continue
                delta = state_journal.finish()
                commit_locus_delta(delta, ref_features_dict, tree_dict)
                if delta.score_text:
                    fw_score.write(delta.score_text)
                step8_genes_emitted += 1
                if getattr(args, "miniprot_rescue", False):
                    # Iter 23: record default Step-8 ref genes for the rescue dedup.
                    emitted_ref_gene_ids.add(lifton_gene.ref_gene_id)
            except Exception as e:
                error_text = safe_exception_text(e)
                logger.log_error(
                    "Error during miniprot text output serialization "
                    f"({mtrans.id}): {error_text}"
                )
                _record_pipeline_failure(
                    args, "process_miniprot_loci", error_text,
                    details={"mRNA": getattr(mtrans, "id", "<unknown>")},
                )

        step8_max_inflight_observed = (
            1 if step8_candidates_processed else 0
        )
        step8_analysis_submitted = step8_candidates_processed

    manifest.record_count(
        "miniprot_candidates_processed", step8_candidates_processed,
    )
    manifest.record_count(
        "miniprot_genes_emitted", step8_genes_emitted,
    )
    manifest.record_count(
        "step8_max_inflight_observed", step8_max_inflight_observed,
    )
    manifest.record_count(
        "miniprot_candidates_pre_overlap_elided",
        step8_pre_overlap_elided,
    )
    manifest.record_count(
        "miniprot_candidates_submitted_for_analysis",
        step8_analysis_submitted,
    )
    manifest.record_count(
        "step8_child_batch_calls", step8_child_batch_calls,
    )
    manifest.record_count(
        "step8_child_batch_fallbacks", step8_child_batch_fallbacks,
    )
    manifest.record_count(
        "step8_child_scalar_materializations",
        step8_child_scalar_materializations,
    )

    # Iteration 23: clean separate-pass miniprot-only rescue. Runs AFTER the
    # Step-8 loop has fully closed, so no rescue mutates tree_dict during a
    # default decision -> the default Step-7+8 output is byte-identical OFF-vs-ON
    # and ON = default ∪ rescues ⊇ default = OFF (off ⊆ on by construction).
    # Flag-gated -> the default path never imports the module.
    if getattr(args, "miniprot_rescue", False):
        from lifton import miniprot_rescue as _miniprot_rescue
        args._rescue_failure_records = []
        rescued_genes = _miniprot_rescue.rescue_miniprot_only_pass(
            m_feature_db, ref_db, tree_dict, tgt_fai, ref_proteins, ref_trans,
            ref_features_dict, m_id_2_ref_id_trans_dict, ref_features_len_dict,
            ref_trans_exon_num_dict, ref_features_reverse_dict,
            emitted_ref_gene_ids, fw, fw_score, transcripts_stats_dict, args)
        manifest.record_count("miniprot_rescued_genes", rescued_genes)
        for failure in args._rescue_failure_records:
            _record_pipeline_failure(
                args, "miniprot_rescue",
                failure.get("message", "miniprot rescue failed"),
                details={
                    key: value for key, value in failure.items()
                    if key != "message"
                },
            )

    t12 = time.process_time()
    ################################
    # Step 9: Printing stats
    ################################
    _switch_manifest_phase(args, "write_reports")
    try:
        stats.print_report(ref_features_dict, transcripts_stats_dict, fw_unmapped, fw_extra_copy, fw_mapped_feature, fw_mapped_trans, debug=args.debug, fw_feature_type=fw_feature_type)
    except Exception as e:
        logger.log_error(f"Failed to print report: {e}")
        _record_pipeline_failure(
            args, "write_reports", e, fatal=True,
        )
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
    # Phase D: cross-locus miniprot rescue (opt-in, post-write rewrite)
    ################################
    # Runs AFTER the writer + score.txt close, so the default Step-7/8 output is
    # complete + untouched when OFF. Gated -> the default path never imports the
    # module. Replaces a WEAK DNA-lifted gene with a strictly-better miniprot
    # model on a DIFFERENT chromosome (drop-the-block-then-append => no dup ids).
    if getattr(args, 'cross_locus_rescue', False) or \
            os.environ.get("LIFTON_CROSS_LOCUS_RESCUE", "0").strip().lower() \
            not in ("", "0", "false", "no", "off"):
        _switch_manifest_phase(args, "cross_locus_rescue")
        from lifton import cross_locus_rescue as _xlocus
        _xlocus.cross_locus_rescue_pass(
            staged_output_path, f"{lifton_outdir}/score.txt", m_feature_db, ref_db,
            tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict,
            m_id_2_ref_id_trans_dict, ref_features_len_dict,
            ref_trans_exon_num_dict, ref_features_reverse_dict, args)

    ################################
    # Step 10: Validate output GFF3
    ################################
    _switch_manifest_phase(args, "validate_output_structure")
    structural_result = gff3_validator.validate_gff3_structure(
        staged_output_path,
    )
    manifest.record_validation(
        "gff3_structural",
        passed=structural_result.is_valid,
        errors=len(structural_result.errors),
        warnings=len(structural_result.warnings),
        details={"path": staged_output_path},
    )
    if not structural_result.is_valid:
        gff3_validator.print_validation_report(
            structural_result,
            verbose=True,
        )
        error = LiftOnValidationError(
            "output GFF3 failed mandatory structural validation with "
            f"{len(structural_result.errors)} error(s)"
        )
        _record_pipeline_failure(
            args,
            "validate_output_structure",
            error,
            fatal=True,
            details={"errors": len(structural_result.errors)},
        )
        partial_path = output_transaction.abort()
        logger.log_error(
            f"Invalid output was preserved at {partial_path}"
        )
        _finish_run_manifest(args, "failed")
        raise error

    if getattr(args, 'validate_output', False):
        _switch_manifest_phase(args, "validate_output")
        print("\n\n*********************************************", file=sys.stderr)
        print("** Validating output GFF3                  **", file=sys.stderr)
        print("*********************************************", file=sys.stderr)
        val_result = gff3_validator.validate_gff3_file(
            gff3_path=staged_output_path,
            check_hierarchy=True,
            check_cds_phase=True,
            check_containment=True,
            check_lifton_attrs=True,
        )
        verbose = getattr(args, 'validate_verbose', False)
        gff3_validator.print_validation_report(val_result, verbose=verbose)
        manifest.record_validation(
            "gff3",
            passed=val_result.is_valid,
            errors=len(val_result.errors),
            warnings=len(val_result.warnings),
            details={"path": staged_output_path},
        )
        if not val_result.is_valid:
            print(
                f"\n[LiftOn] Output GFF3 has {len(val_result.errors)} error(s). "
                f"See report above for details.",
                file=sys.stderr,
            )
            error = LiftOnValidationError(
                f"output GFF3 failed validation with "
                f"{len(val_result.errors)} error(s)"
            )
            _record_pipeline_failure(
                args, "validate_output", error, fatal=True,
                details={"errors": len(val_result.errors)},
            )
            partial_path = output_transaction.abort()
            logger.log_error(
                f"Invalid output was preserved at {partial_path}"
            )
            _finish_run_manifest(args, "failed")
            raise error

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
        console = _console_stream(args)
        print("Time taken for each step:", file=console)
        print(f"Reading target & reference genomes: {reading_target_reference_genomes}", file=console)
        print(f"Creating reference annotation database: {creating_reference_annotation_database}", file=console)
        print(f"Get all reference features to liftover: {get_all_reference_features_to_lift}", file=console)
        print(f"Extract protein & DNA dictionaries: {extract_protein_dna_dictionaries}", file=console)
        print(f"Run liftoff & miniprot: {run_liftoff_miniprot}", file=console)
        print(f"Create liftoff database: {create_liftoff_database}", file=console)
        print(f"Create miniprot database: {create_miniprot_database}", file=console)
        print(f"Miniprot 2 Liftoff ID mapping: {miniprot_2_liftoff_id_mapping}", file=console)
        print(f"Process Liftoff genes & transcripts: {process_liftoff_genes_transcripts}", file=console)
        print(f"Process miniprot transcripts: {process_miniprot_transcripts}", file=console)
        print(f"Report stats: {report_stats}", file=console)
        print(f"Overall time: {overall_time}", file=console)
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

    _switch_manifest_phase(args, "publish_output")
    reference_feature_ids = [
        feature_id for feature_id in ref_features_dict
        if feature_id != "LiftOn-gene"
    ]
    mapped_reference_features = sum(
        1 for feature_id in reference_feature_ids
        if ref_features_dict[feature_id].copy_num > 0
    )
    emitted_feature_copies = sum(
        max(0, ref_features_dict[feature_id].copy_num)
        for feature_id in reference_feature_ids
    )
    emitted_transcripts = sum(
        count for category in transcripts_stats_dict.values()
        for count in category.values()
    )
    manifest.record_count("reference_features", len(reference_feature_ids))
    manifest.record_count(
        "mapped_reference_features", mapped_reference_features,
    )
    manifest.record_count(
        "emitted_feature_copies", emitted_feature_copies,
    )
    manifest.record_count(
        "emitted_transcripts", emitted_transcripts,
    )
    manifest.record_count("processed_features", processed_features)
    if reference_feature_ids and emitted_feature_copies == 0:
        _record_pipeline_failure(
            args, "publish_output",
            "no feature hierarchies were emitted from a non-empty reference",
            details={"reference_features": len(reference_feature_ids)},
        )
    manifest.record_count(
        "pipeline_failures", len(getattr(args, "_pipeline_failures", []))
    )
    failures = getattr(args, "_pipeline_failures", [])
    if failures and not getattr(args, "allow_partial_output", False):
        error = LiftOnPartialOutputError(
            f"{len(failures)} locus/report failure(s) prevented publication; "
            "rerun with --allow-partial-output to publish the staged result"
        )
        _record_pipeline_failure(
            args, "publish_output", error, fatal=True,
            details={"failure_count": len(failures)},
            affects_completeness=False,
        )
        partial_path = output_transaction.abort()
        logger.log_error(f"Partial output was preserved at {partial_path}")
        _finish_run_manifest(args, "failed")
        raise error

    output_transaction.commit()
    manifest.set_backend_choice(
        "output",
        "stdout" if args.output == "stdout" else os.path.abspath(args.output),
    )
    _switch_manifest_phase(args, None)
    _finish_run_manifest(
        args, "partial_success" if failures else "success"
    )


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
    if all(hasattr(args, name) for name in (
            "output", "target", "reference", "reference_annotation")):
        outdir, lifton_outdir = _run_outdirs(args.output)
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(lifton_outdir, exist_ok=True)
        _ensure_run_manifest(args, outdir, lifton_outdir)
    try:
        # Only preflight tools the selected execution path will actually launch.
        # Evaluation and valid precomputed -L/-M inputs are intentionally usable
        # on machines without the aligner binaries.
        if not getattr(args, "evaluation", False):
            needs_miniprot = not (
                getattr(args, "miniprot", None) is not None
                and os.path.exists(args.miniprot)
            )
            if (needs_miniprot
                    and not run_miniprot.check_miniprot_installed()):
                sys.exit(
                    "miniprot is not installed. Please install miniprot before "
                    "running LiftOn, or provide a valid precomputed -M file."
                )

            needs_minimap2 = not (
                getattr(args, "liftoff", None) is not None
                and os.path.exists(args.liftoff)
            )
            if needs_minimap2:
                native_requested = (
                    bool(getattr(args, "native", False))
                    and bool(os.environ.get("LIFTON_NATIVE_LIFTOFF_ALIGN"))
                )
                if native_requested:
                    from lifton.native_bindings import is_mappy_available
                    needs_minimap2 = not is_mappy_available()
            custom_minimap2 = getattr(args, "m", None)
            minimap2_path = custom_minimap2 or "minimap2"
            minimap2_ok = True
            if needs_minimap2:
                minimap2_ok = (
                    run_liftoff.check_minimap2_installed(custom_minimap2)
                    if custom_minimap2
                    else run_liftoff.check_minimap2_installed()
                )
            if needs_minimap2 and not minimap2_ok:
                sys.exit(
                    f"minimap2 is not installed or runnable at "
                    f"'{minimap2_path}'. Install it, pass -m PATH, provide a "
                    "valid precomputed -L file, or use the supported native "
                    "Liftoff path."
                )
        if getattr(args, "output", None) == "stdout":
            args._gff_stdout = sys.stdout
            with redirect_stdout(sys.stderr):
                run_all_lifton_steps(args)
        else:
            run_all_lifton_steps(args)
    except BaseException as exc:
        # Preserve staged GFF3 bytes for diagnosis, but never replace a valid
        # destination after an interrupted/failed run.
        transaction = getattr(args, "_output_transaction", None)
        if transaction is not None:
            try:
                partial_path = transaction.abort()
                logger.log_error(
                    f"Unpublished output was preserved at {partial_path}"
                )
            except Exception:
                # Already committed, or the staging filesystem itself failed.
                pass

        manifest = getattr(args, "_run_manifest", None)
        if manifest is not None:
            status = manifest.to_dict()["run"]["status"]
            if status == "running":
                manifest.record_failure("pipeline", exc, fatal=True)
                _finish_run_manifest(args, "failed")
            else:
                manifest.write(args._run_manifest_path)
        if isinstance(exc, (LiftOnPartialOutputError,
                            LiftOnValidationError)):
            raise SystemExit(2) from None
        raise
