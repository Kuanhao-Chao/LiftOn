import re, sys, os, copy
from Bio.Seq import Seq
from lifton import align, lifton_class, run_liftoff, run_miniprot, logger
from lifton.liftoff import liftoff_main
# Iteration 16: the three pure helpers below moved to the dependency-free
# leaf module lifton.coreutils to break the lifton_utils <-> lifton_class
# import cycle. Re-exported here so lifton_utils.<helper> keeps resolving.
from lifton import coreutils
from lifton.tool_execution import EXTERNAL_ALIGNER_INSTALL_HELP
import gffutils
from lifton.coreutils import (  # noqa: F401
    custom_bisect_insert,
    get_ID_base,
    segments_overlap_length,
)

def check_miniprot_installed():
    """
        This function checks if miniprot is installed.

        Parameters:
        None

        Returns:
        Exist the program if miniprot is not installed.
    """
    miniprot_installed = run_miniprot.check_miniprot_installed()
    # print("miniprot_installed: ", miniprot_installed)
    if not miniprot_installed:
        sys.exit("miniprot is not installed. " + EXTERNAL_ALIGNER_INSTALL_HELP)


def get_truncated_protein(ref_proteins):
    """
        This function gets the truncated proteins.

        Parameters:
        - ref_proteins: reference proteins dictionary

        Returns:
        truncated_proteins: truncated proteins dictionary
    """
    truncated_proteins = {}
    for record in ref_proteins.keys():
        protein = ref_proteins[record]
        if not check_protein_valid(str(protein)):
            truncated_proteins[record] = protein
    # print("truncated_proteins: ", len(truncated_proteins))
    # print("good_protein: ", good_protein)
    # print("bad_protein: ", bad_protein)
    return truncated_proteins


def count_truncated_proteins(ref_proteins):
    """Number of reference proteins that fail ``check_protein_valid``.

    The pipeline only ever reported ``len()`` of the dict :func:`
    get_truncated_protein` builds, so at human scale it materialised a
    dictionary of pyfaidx records -- and kept it alive for the whole run --
    to produce a single integer.
    """
    return sum(
        1 for record in ref_proteins.keys()
        if not check_protein_valid(str(ref_proteins[record]))
    )


def write_seq_2_file(outdir, ref_seqs, target):
    """
        This function writes the reference sequences to a file.

        Parameters:
        - outdir: output directory
        - ref_seqs: reference sequences dictionary
        - target: target type ('truncated_proteins', 'proteins', or 'transcripts')

        Returns:
        ref_seqs_file: reference sequences file path
    """
    if target == "truncated_proteins":
        ref_seqs_file = outdir + "/proteins_truncated.fa"
    elif target == "proteins":
        ref_seqs_file = outdir + "/proteins.fa"
    elif target == "transcripts":
        ref_seqs_file = outdir + "/transcripts.fa"
    fw = open(ref_seqs_file, 'w')
    # Iterate through the original FASTA and write the records to the new FASTA file
    for record in ref_seqs.keys():
        seq = ref_seqs[record]
        fw.write(f'>{record}\n{seq}\n')
    fw.close()
    return ref_seqs_file


def check_protein_valid(protein):
    """
        This function checks if the protein is valid.

        Parameters:
        - protein: protein sequence

        Returns:
        True if the protein is valid, False otherwise.
    """
    # The length of the protein has to be greater than 0
    if len(protein) == 0:
        return False
    # Start with M
    if protein[0] != "M":
        return False
    # End with *
    if protein[-1] != "*":
        return False
    # Only 1 * in the string
    if protein.count("*") != 1:
        return False
    return True
    

def exec_liftoff(outdir, ref_db, args):
    """
        This function executes liftoff.

        Parameters:
        - outdir: output directory
        - args: arguments

        Returns:
        liftoff_annotation: liftoff annotation file path
    """
    liftoff_annotation = args.liftoff
    if liftoff_annotation is None or not os.path.exists(liftoff_annotation):
        print("\n*********************")
        print("** Running Liftoff **")
        print("*********************")
        liftoff_annotation = run_liftoff.run_liftoff(outdir, ref_db, args)
    return liftoff_annotation


def exec_miniprot(outdir, args, tgt_genome, ref_proteins_file):
    """
        This function executes miniprot.

        Parameters:
        - outdir: output directory
        - args: arguments
        - tgt_genome: target genome
        - ref_proteins_file: reference protein file

        Returns:
        miniprot_annotation: miniprot annotation file path
    """
    miniprot_annotation = args.miniprot
    if miniprot_annotation is None or not os.path.exists(miniprot_annotation):
        # Preflight ONLY when miniprot is actually going to run. Checking before the
        # -M short-circuit made the documented cached workflow
        #   lifton -g ref.gff3 -L prior.gff3 -M prior_miniprot.gff3 ...
        # exit 1 on a machine without miniprot installed, even though the binary is
        # never invoked. (main() already gates its own preflight the same way.)
        check_miniprot_installed()
        print("\n**********************")
        print("** Running miniprot **")
        print("**********************")
        miniprot_annotation = run_miniprot.run_miniprot(outdir, args, tgt_genome, ref_proteins_file)
    return miniprot_annotation


def get_ID(feature):
    """
        This function gets the ID.

        Parameters:
        - feature: feature

        Returns:    
        the original ID and the ID base
    """
    # print("feature: ", feature)
    # id_spec={"ID"}
    id = feature.id
    id_base = get_ID_base(id)
    return id, id_base


def get_parent_features_to_lift(feature_types_file):
    """
        This function gets the parent features to lift.

        Parameters:
        - feature_types_file: feature types file (from '-f' / '--features')

        Returns:
        feature_types: feature types
    """
    feature_types = ["gene"]
    if feature_types_file is not None:
        feature_types = []
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types


def _flat_fallback_types():
    """Top-level types worth lifting when no type bears a hierarchy at all.

    A positive list, not an exclusion list: a flat annotation's genes are its
    `CDS` rows (bakta output, a miniprot GFF), while its other top-level rows
    are landmarks and regulatory features that were never lift targets. An
    exclusion list would have to anticipate every one of those -- and a first
    attempt at one selected childless `enhancer` rows, which an existing test
    caught.
    """
    from lifton.gff3_validator import GENE_TYPES, TRANSCRIPT_TYPES
    return frozenset(GENE_TYPES) | frozenset(TRANSCRIPT_TYPES) | {"CDS"}


def get_gene_like_feature_types(ref_db, sample_cap=5000):
    """Auto-detect the "gene-like" top-level parent feature types in the
    reference annotation: every feature type that has at least one TOP-LEVEL
    instance (no Parent) bearing a child feature (a transcript/exon hierarchy).

    This generalises the lift beyond the hardcoded `["gene"]` default to
    pseudogenes, ncRNA_gene, structured mobile elements, etc., in an
    annotation-source-agnostic way (Iteration 5 `--lift-gene-like`). Childless
    meta/regulatory features (region, enhancer, promoter, match, …) and pure
    child types (CDS, exon, mRNA, …, which carry a Parent) are excluded.

    `sample_cap` bounds the per-type scan so a child-heavy type (e.g. CDS, whose
    instances all carry Parent and are skipped) cannot trigger a full-DB walk;
    realistic gene-like types surface a child-bearing top-level instance far
    within the cap.

    When nothing qualifies, fall back to the top-level gene, transcript or `CDS`
    types the annotation actually has. A prokaryotic or otherwise flat
    annotation -- bakta output, a miniprot GFF -- has top-level `CDS` rows and
    no `gene` at all, and the old `["gene"]` fallback selected nothing from it,
    so the run died several steps later inside vendored Liftoff (GH #37). Falls
    back to `["gene"]` when the annotation has no such type either.

    The fallback list applies ONLY when nothing bears a hierarchy. Filtering the
    detection itself would drop a type that genuinely has one -- `match` over
    `match_part`, say -- from annotations the old code lifted.
    """
    conn = ref_db.db_connection
    liftable = _flat_fallback_types()
    gene_like, top_level = set(), set()
    for ftype in conn.featuretypes():
        for i, locus in enumerate(conn.features_of_type(ftype)):
            if i >= sample_cap:
                break
            if "Parent" in locus.attributes:
                continue                      # child-level instance; keep looking
            if ftype in liftable:
                top_level.add(ftype)
            if any(True for _ in conn.children(locus, level=1)):
                gene_like.add(ftype)
                break                         # one child-bearing top-level is enough
    if gene_like:
        return sorted(gene_like)
    return sorted(top_level) if top_level else ["gene"]


def LiftOn_eval_alignment(eval_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status):
    eval_aln = align.lifton_parasail_align(eval_trans, locus, tgt_fai, ref_proteins, ref_trans_id)
    if eval_aln != None:
        lifton_status.lifton_aa = eval_aln.identity
    return eval_aln


def LiftOn_liftoff_alignment(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status):
    liftoff_aln = align.lifton_parasail_align(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id)
    if liftoff_aln != None:
        lifton_status.liftoff = liftoff_aln.identity
    return liftoff_aln


def LiftOn_miniprot_alignment(chromosome, transcript, m_id_dict, m_feature_db, tree_dict, fai, ref_proteins, ref_trans_id, lifton_status):
    """
        This function checks the miniprot alignment.

        Parameters:
        - chromosome: chromosome
        - transcript: transcript gffutils feature
        - lifton_status: Lifton_Status instance
        - m_id_dict: miniprot id dictionary
        - m_feature_db: miniprot feature database
        - tree_dict: tree dictionary
        - fai: reference fasta index
        - ref_proteins: reference proteins dictionary
        - ref_trans_id: reference transcript ID

        Returns:
        m_lifton_aln: miniprot lifton alignment
        has_valid_miniprot: True if the miniprot transcript is valid, False otherwise.
    """
    m_lifton_aln = None
    has_valid_miniprot = False
    if (ref_trans_id in m_id_dict.keys()) and (ref_trans_id in ref_proteins.keys()):
        m_ids = m_id_dict[ref_trans_id]
        for m_id in m_ids:
            ##################################################
            # Check 1: Check if the miniprot transcript is overlapping the current gene locus
            ##################################################
            m_entry = m_feature_db[m_id]
            _, overlap = segments_overlap_length((m_entry.start, m_entry.end), (transcript.start, transcript.end))
            if not overlap or m_entry.seqid != transcript.seqid:
                # "Not overlapped"
                continue
            ##################################################
            # Check 2: reference overlapping status
            #   1. Check it the transcript overlapping with the next gene
            # Check the miniprot protein overlapping status
            # The case I should not process the transcript 
            #   1. The Liftoff does not overlap with other gene
            #   2. The miniprot protein overlap the other gene
            ##################################################
            ovps_liftoff = tree_dict[chromosome].overlap(transcript.start, transcript.end)
            ovps_miniprot = tree_dict[chromosome].overlap(m_entry.start, m_entry.end)
            miniprot_cross_gene_loci = False
            liftoff_set = set()
            for ovp_liftoff in ovps_liftoff:
                liftoff_set.add(ovp_liftoff[2])
            for ovp_miniprot in ovps_miniprot:
                if ovp_miniprot[2] not in liftoff_set:
                    # Miniprot overlap to more genes
                    miniprot_cross_gene_loci = True
                    break
            if miniprot_cross_gene_loci:
                continue
            # Valid miniprot transcript exists => check if the miniprot transcript is valid
            has_valid_miniprot = True
            miniprot_trans = lifton_class.Lifton_TRANS(m_id, "", "", 0, m_entry, {})
            # ONE query, not two identical ones: this exact statement used to
            # run twice per candidate, 7,038 times on drosophila -- 9.6% of every
            # SQL statement Step 7 issued. The second consumer needs its own
            # objects because `Lifton_EXON.__init__` rewrites `featuretype`, so
            # feeding it the exons' features would hand `add_cds` rows already
            # relabelled 'exon'. Cloning costs ~5 us against a SQL round-trip.
            children = coreutils.drop_redundant_stop_codons(
                m_feature_db.children(
                    m_entry, featuretype=('CDS', 'stop_codon'),
                    order_by='start'))
            # miniprot repeats the terminal CDS's stop codon as a nested
            # `stop_codon` row. Every one of those rows is already covered by a
            # sibling CDS, so ingesting it added a 3 bp exon inside the terminal
            # exon and let `add_cds` overwrite the real terminal CDS. A stop
            # codon that genuinely sits outside the CDS is still kept.
            for exon in children:
                miniprot_trans.add_exon(exon)
            cds_num = 0
            for cds in children:
                cds_num += 1
                miniprot_trans.add_cds(coreutils.clone_feature(cds))
            tmp_m_lifton_aln = align.lifton_parasail_align(miniprot_trans, m_entry, fai, ref_proteins, ref_trans_id)
            if tmp_m_lifton_aln is None:
                # `lifton_parasail_align` returns None when it cannot build a protein
                # alignment -- e.g. a miniprot mRNA with no CDS/stop_codon children, so
                # the coding sequence (and therefore `protein_seq`) is empty. The old
                # code assigned it and then read `.identity`, raising AttributeError;
                # the caller catches that per-locus, so ONE unusable miniprot candidate
                # silently dropped the ENTIRE Liftoff gene from the output even though
                # its DNA lift was fine. Skip the candidate instead.
                continue
            if m_lifton_aln is None or tmp_m_lifton_aln.identity > lifton_status.miniprot:
                m_lifton_aln = tmp_m_lifton_aln
                lifton_status.miniprot = m_lifton_aln.identity
    # Only advertise a valid miniprot alignment when one was actually produced: callers
    # gate `chaining_algorithm(...)` on this flag and would dereference a None alignment.
    return m_lifton_aln, (has_valid_miniprot and m_lifton_aln is not None)


def _parent_type(ref_db, parent_id):
    """Iteration 20: return the featuretype of the feature `parent_id`, or None
    if it cannot be resolved. Used by get_ref_liffover_features to skip a child
    instance whose parent is itself a lifted type (see extract_sequence's
    parent_is_listed_type for the full rationale)."""
    try:
        return ref_db.db_connection[parent_id].featuretype
    except Exception:
        return None


def get_ref_liffover_features(features, ref_db, intermediate_dir, args):
    """
        This function gets the reference liftover features.

        Parameters:
        - features: list of features to liftover
        - ref_db: reference database

        Returns:
        ref_features_dict: reference features dictionary (gene id -> transcript id)
        ref_features_reverse_dict: reference features reverse dictionary (transcript id -> gene id)
    """

    fw_gene = open(f"{intermediate_dir}/ref_feature.txt", "w")
    fw_trans = open(f'{intermediate_dir}/ref_transcript.txt', 'w')    
    ref_features_dict = {}
    ref_features_len_dict = {}
    ref_features_reverse_dict = {}
    ref_trans_exon_num_dict = {}
    new_gene_feature = lifton_class.Lifton_feature("Lifton-gene")
    ref_features_dict["LiftOn-gene"] = new_gene_feature
    feature_set = set(features)
    for f_itr in features:
        for locus in ref_db.db_connection.features_of_type(f_itr):
            # Iteration 20: skip a child instance whose parent is itself one of
            # the lifted types — it is lifted under its parent gene (the level-1
            # transcripts loop below), so adding it here as its own locus would
            # emit a duplicate gene model. Mirrors extract_sequence's
            # parent_is_listed_type so extraction and lift stay consistent;
            # inlined here to avoid a new import edge into lifton_utils. No-op on
            # top-level-only annotations (parent type not in the lift set).
            parents = locus.attributes.get("Parent")
            if parents and any(
                _parent_type(ref_db, pid) in feature_set for pid in parents):
                continue
            # ONE recursive CDS query per locus, answering both consumers: the
            # coding test below needs only a count, and ref_features_len_dict
            # needs the ordered extremes. This used to be the same recursive
            # query issued twice -- once unordered here and once ordered below
            # -- which is the Iteration-18 pattern collapsed in Step 3 and
            # never applied here. Ordering cannot change a count, so the
            # coding test sees exactly what it saw before.
            CDS_children = list(ref_db.db_connection.children(
                locus, featuretype='CDS', order_by='start'))
            feature = lifton_class.Lifton_feature(locus.id)
            feature.feature_type = locus.featuretype
            # Write out reference gene features IDs
            # Decide if its type
            # Check for gene type/biotype attribute in order of preference
            gene_type_key = None
            if args.annotation_database.upper() == "REFSEQ":
                # For RefSeq, check gene_biotype first, then biotype as fallback
                if "gene_biotype" in locus.attributes.keys():
                    gene_type_key = "gene_biotype"
                elif "biotype" in locus.attributes.keys():
                    gene_type_key = "biotype"
            elif args.annotation_database.upper() == "GENCODE" or args.annotation_database.upper() == "ENSEMBL" or args.annotation_database.upper() == "CHESS":
                # For GENCODE/ENSEMBL/CHESS, check gene_type first, then biotype as fallback
                if "gene_type" in locus.attributes.keys():
                    gene_type_key = "gene_type"
                elif "gene_biotype" in locus.attributes.keys():
                    gene_type_key = "gene_biotype"
                elif "biotype" in locus.attributes.keys():
                    gene_type_key = "biotype"
            else:
                # For other databases, try all possible keys in order
                if "gene_biotype" in locus.attributes.keys():
                    gene_type_key = "gene_biotype"
                elif "gene_type" in locus.attributes.keys():
                    gene_type_key = "gene_type"
                elif "biotype" in locus.attributes.keys():
                    gene_type_key = "biotype"

            if gene_type_key is not None:
                feature.biotype = locus.attributes[gene_type_key][0]
                if locus.attributes[gene_type_key][0] == "protein_coding" and len(CDS_children) > 0:
                    feature.is_protein_coding = True
                    fw_gene.write(f"{locus.id}\tcoding\n")
                elif (locus.attributes[gene_type_key][0] == "lncRNA" or locus.attributes[gene_type_key][0] == "ncRNA"):
                    feature.is_non_coding = True
                    fw_gene.write(f"{locus.id}\tnon-coding\n")
                else:
                    fw_gene.write(f"{locus.id}\tother\n")
            else:
                fw_gene.write(f"{locus.id}\tother\n")
            exon_children = list(ref_db.db_connection.children(locus, featuretype='exon', level=1, order_by='start'))
            if len(exon_children) > 0:
                __process_ref_liffover_features(locus, ref_db, None)
            else:
                transcripts = ref_db.db_connection.children(locus, level=1)
                for transcript in list(transcripts):
                    __process_ref_liffover_features(transcript, ref_db, feature)
                    ref_features_reverse_dict[transcript.id if not args.evaluation_liftoff_chm13 else locus.id[4:]] = locus.id
                    all_CDS_in_trans = list(ref_db.db_connection.children(transcript, featuretype='CDS', order_by='start'))
                    if len(all_CDS_in_trans) > 0:
                        ref_trans_exon_num_dict[transcript.id if not args.evaluation_liftoff_chm13 else locus.id[4:]] = len(all_CDS_in_trans)
                    else:
                        ref_trans_exon_num_dict[transcript.id if not args.evaluation_liftoff_chm13 else locus.id[4:]] = 0
                    # Write out reference trans feature IDs
                    if feature.is_protein_coding and transcript.featuretype == "mRNA":
                        fw_trans.write(f"{transcript.id}\tcoding\n")
                    elif feature.is_non_coding and (transcript.featuretype == "ncRNA" or transcript.featuretype == "nc_RNA" or transcript.featuretype == "lncRNA" or transcript.featuretype == "lnc_RNA"):
                        fw_trans.write(f"{transcript.id}\tnon-coding\n")
                    else:
                        fw_trans.write(f"{transcript.id}\tother\n")
            ref_features_dict[locus.id if not args.evaluation_liftoff_chm13 else locus.id[5:]] = feature
            all_CDS_children = CDS_children
            if len(all_CDS_children) > 0:
                ref_features_len_dict[locus.id if not args.evaluation_liftoff_chm13 else locus.id[5:]] = all_CDS_children[-1].end - all_CDS_children[0].start + 1
            else:
                ref_features_len_dict[locus.id if not args.evaluation_liftoff_chm13 else locus.id[5:]] = 0
    fw_gene.close()
    fw_trans.close()
    return ref_features_dict, ref_features_len_dict, ref_features_reverse_dict, ref_trans_exon_num_dict


def __process_ref_liffover_features(locus, ref_db, feature):
    if feature != None:
        feature.children.add(locus.id)


def miniprot_id_mapping(m_feature_db):
    """
        This function creates a dictionary of miniprot id to reference id.

        Parameters:
        - m_feature_db: miniprot feature database

        Returns:
        ref_id_2_m_id_trans_dict: reference id to miniprot transcript ids dictionary
        m_id_2_ref_id_trans_dict: miniprot transcript id to reference id dictionary
    """
    ref_id_2_m_id_trans_dict = {}
    m_id_2_ref_id_trans_dict = {}
    if m_feature_db is None:
        return ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict
        
    for feature in m_feature_db.features_of_type("mRNA"):
        miniprot_id = feature["ID"][0]
        aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]
        if aa_trans_id in ref_id_2_m_id_trans_dict.keys():
            ref_id_2_m_id_trans_dict[aa_trans_id].append(miniprot_id)
        else:
            ref_id_2_m_id_trans_dict[aa_trans_id] = [miniprot_id]
        m_id_2_ref_id_trans_dict[miniprot_id] = aa_trans_id
    return ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict


def get_ref_ids_liftoff(ref_features_dict, liftoff_gene_id, liftoff_trans_id):
    """
        This function gets the reference IDs from Liftoff IDs.

        Parameters:
        - ref_features_dict: reference features dictionary
        - liftoff_gene_id: Liftoff gene ID
        - liftoff_trans_id: Liftoff transcript ID

        Returns:
        # Three cases that will be passed into this function
        1.  Liftoff_gene_id, None
        2.  None, Liftoff_trans_id
        3.  Liftoff_gene_id, Liftoff_trans_id
    """
    if liftoff_gene_id is None:
        ref_trans_id = __extract_ref_ids(ref_features_dict, liftoff_trans_id)
        return None, ref_trans_id 

    if liftoff_trans_id is None:
        ref_gene_id = __extract_ref_ids(ref_features_dict, liftoff_gene_id)
        return ref_gene_id, None

    ref_gene_id = __extract_ref_ids(ref_features_dict, liftoff_gene_id)
    if ref_gene_id is None:
        return None, None
    else:
        # ``ref_features_dict`` is keyed by GENE id, so it cannot confirm a
        # transcript-level copy suffix, and passing it here would reject every
        # transcript. Passing ``None`` instead used to disable the check
        # altogether -- ``get_ID_base`` then returned the id verbatim, so a
        # Liftoff ``-copies`` transcript (``rna-X_1``) never resolved to its
        # reference (``rna-X``) and the caller dropped the transcript plus all
        # of its exons and CDS, leaving a bare gene line in the output.
        #
        # The suffix is therefore resolved against the reference annotation
        # itself, by ``resolve_ref_trans_id`` at the call site: it is the only
        # source that can distinguish a copy suffix from an id that genuinely
        # ends in ``_<int>``. This function returns the id as Liftoff wrote it
        # and leaves that resolution to the caller.
        return ref_gene_id, liftoff_trans_id


def __extract_ref_ids(ref_features_dict, liftoff_id):
    if liftoff_id in ref_features_dict.keys():
        return liftoff_id
    else:
        # Pass ref_features_dict to get_ID_base so it can verify the base ID exists
        ref_id = get_ID_base(liftoff_id, ref_features_dict)
        if ref_id in ref_features_dict.keys():
            return ref_id
        else:
            return None


def copy_suffix_base(feature_id):
    """Return the candidate base of a Liftoff ``-copies`` id, or ``None``.

    Liftoff appends ``_<extra_copy_number>`` to the ID (and Parent) of every
    feature of an extra gene copy -- ``liftoff/write_new_gff.py:edit_copy_ids``.
    The base returned here is only ever a *candidate*: reference ids that
    genuinely end in ``_<int>`` exist (``...mrna.FMUND_1``), so callers must try
    the exact id first and fall back to this only on a miss.
    """
    if not feature_id:
        return None
    base, separator, suffix = str(feature_id).rpartition("_")
    if separator and base and suffix.isdigit():
        return base
    return None


def resolve_ref_trans_id(ref_db, candidate, ref_gene_id=None):
    """Resolve a Liftoff transcript id to its reference transcript id.

    Exact id first, then -- and only on a miss -- the Liftoff ``-copies``
    ``_<N>`` base. Trying the exact id first is what makes this safe for
    reference ids that legitimately end in ``_<int>``: such an id is present in
    the reference, so it matches exactly and is never stripped. This mirrors the
    rule ``run_evaluation._copy_base`` already applies on the evaluation path.

    When ``ref_gene_id`` is known, a stripped base is accepted only if it is a
    child of that reference gene, so a copy suffix can never resolve onto some
    other gene's transcript. ``ref_gene_id`` is ``None`` for the 3-level
    hierarchy case (gene -> primary_transcript -> miRNA), where the caller
    resolves from the transcript id alone; the guard is skipped there.

    Returns the resolved reference transcript id, or ``None`` if neither the
    exact id nor a verified copy base exists in the reference annotation.
    """
    # No early-out on a None candidate. The pre-existing behaviour was a plain
    # `ref_db[ref_trans_id]` probe, which this replaces, and the caller can still
    # reach here with None (the reference gene did not resolve). Short-circuiting
    # would stop touching `ref_db` at all on that path, and
    # `tests/test_vulnerabilities.py::TestV1_1a_BareExceptInRefDbLookup` pins
    # that an interrupt raised *by the lookup* still propagates. The probe is
    # harmless -- every backend raises a caught miss for a None key.
    if __ref_db_has(ref_db, candidate):
        return candidate
    base = copy_suffix_base(candidate)
    if base is None or not __ref_db_has(ref_db, base):
        return None
    if ref_gene_id is not None and not __ref_parent_is(ref_db, base, ref_gene_id):
        return None
    return base


# Exactly the pair the pre-existing `ref_db[ref_trans_id]` probe in
# `run_liftoff.process_liftoff` caught, and the same pair `run_evaluation`
# uses. KeyError is the gffbase / dict-style miss and what
# `locus_pipeline._RefDbProxy` raises for an id its parent thread did not
# pre-fetch; FeatureNotFoundError is gffutils. Anything else is a real backend
# failure and must keep propagating -- catching it here would turn a broken
# database into a silent genome-wide drop, which is the shape of the bug this
# whole change exists to fix.
_MISSING_FEATURE_ERRORS = (KeyError, gffutils.exceptions.FeatureNotFoundError)


def __ref_db_has(ref_db, feature_id):
    """True if ``ref_db[feature_id]`` resolves."""
    try:
        ref_db[feature_id]
    except _MISSING_FEATURE_ERRORS:
        return False
    return True


def __ref_parent_is(ref_db, feature_id, ref_gene_id):
    """True if ``feature_id``'s reference Parent includes ``ref_gene_id``."""
    try:
        parents = ref_db[feature_id].attributes.get("Parent") or []
    except _MISSING_FEATURE_ERRORS:
        return False
    if isinstance(parents, str):
        parents = [parents]
    return any(str(parent) == str(ref_gene_id) for parent in parents)


def get_ref_ids_miniprot(ref_features_reverse_dict, miniprot_trans_id, m_id_2_ref_id_trans_dict):
    if miniprot_trans_id not in m_id_2_ref_id_trans_dict.keys():
        return None, None
    ref_trans_id = m_id_2_ref_id_trans_dict[miniprot_trans_id]
    if ref_trans_id not in ref_features_reverse_dict.keys():
        return None, ref_trans_id
    return ref_features_reverse_dict[ref_trans_id], ref_trans_id


def print_lifton_status(transcript_id, transcript, lifton_status, DEBUG=False):
    final_status = ";".join(lifton_status.status)
    logger.log(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton_dna}\t{lifton_status.lifton_aa}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n", debug=DEBUG)


def write_lifton_status(fw_score, transcript_id, transcript, lifton_status):
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton_dna}\t{lifton_status.lifton_aa}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def write_lifton_eval_status(fw_score, transcript_id, transcript, lifton_status):
    """Write the established five-column evaluation record."""
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.lifton_dna}\t{lifton_status.lifton_aa}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def write_lifton_chains(fw_chain, transcript_id, chains):
    chain_ls = ";".join(chains)
    fw_chain.write(f"{transcript_id}\t{chain_ls}\n")


def check_ovps_ratio(mtrans, mtrans_interval, overlap_ratio, tree_dict):
    """
        This function checks the overlap ratio.

        Parameters:
        - mtrans: miniprot transcript gffutils feature
        - mtrans_interval: miniprot transcript interval
        - overlap_ratio: overlap ratio
        - tree_dict: tree dictionary

        Returns:
        True if the overlap ratio is greater than the threshold, False otherwise.
    """
    is_overlapped = False
    if mtrans.seqid not in tree_dict.keys():
        return False
    # Bug fix #3 (Phase 5): IntervalTree.overlap requires (begin, end)
    # ints or an Interval, not a raw tuple. Accept either an
    # intervaltree.Interval or a (begin, end) tuple here.
    if hasattr(mtrans_interval, "begin"):
        begin, end = mtrans_interval.begin, mtrans_interval.end
    else:
        begin, end = mtrans_interval[0], mtrans_interval[1]
    ovps = tree_dict[mtrans.seqid].overlap(begin, end)
    for ovp in ovps:
        ovp_len, _ = segments_overlap_length((mtrans_interval[0], mtrans_interval[1]), (ovp[0], ovp[1]))
        ref_len = ovp[1] - ovp[0] + 1
        target_len = mtrans_interval[1] - mtrans_interval[0] + 1
        # Overlapping does not extend the ratio of the reference
        if (ovp_len / min(ref_len, target_len)) > overlap_ratio:
            is_overlapped = True
            break
    return is_overlapped
