from lifton import align, coreutils, get_id_fraction, variants, logger
from lifton.cds_id_allocator import CdsIdAllocator
from lifton.io import feature_serializer
import copy, os, re
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree


def _containment_normalize_enabled():
    """Iteration-24 GFF3 parent-child containment normalization (default ON).

    LiftOn's chaining + ORF-boundary patching can leave a transcript whose
    child features fall OUTSIDE their parent span (invalid GFF3 the reference
    does not have): a CDS extended a few bp past its exon to capture a stop
    codon on a frameshift-corrected RefSeq model (CDS > exon > mRNA), or a
    rebuilt non-monotonic exon list that inverts the mRNA span (start > end,
    the Iter-21 crash). ``normalize_containment`` (run at write time) extends
    each exon to cover its CDS, sorts the exons, and sets the transcript/gene
    span to the child envelope — guaranteeing valid containment.

    Default ON. Set ``LIFTON_NO_CONTAINMENT_NORMALIZE=1`` to reproduce the
    pre-Iteration-24 bytes (manuscript reproduction / A-B baseline).
    """
    return os.environ.get("LIFTON_NO_CONTAINMENT_NORMALIZE", "").lower() not in ("1", "true", "yes")


# Attributes that IDENTIFY or RE-PARENT a CDS row rather than describe it, and
# so are never carried across a rebuild: ``ID`` belongs to the write-funnel
# allocator (``Lifton_TRANS._normalize_cds_ids``), ``Parent`` is rewritten to the
# emitting transcript, and ``extra_copy_number`` is LiftOn bookkeeping that never
# reaches the output.
_CDS_STRUCTURAL_ATTR_KEYS = ("ID", "Parent", "extra_copy_number")


def _cds_attr_carry_enabled():
    """Carry a CDS's descriptive attributes across a rebuild (default ON).

    The chaining merge and the ORF-rescue boundary patch REBUILD a transcript's
    CDS list, and both used to reset the emitted attributes to ``{Parent}``. The
    result was an internally inconsistent file: on drosophila, 7,127 of 7,142
    pure-Liftoff CDS rows kept ``Dbxref``/``product``/``protein_id``/``gene``/…
    while only 21 of 151,277 merged rows did — even though every merged mRNA row
    kept the very same attributes. Reference CDS segments of one transcript share
    those attributes (a discontinuous feature), so a rebuilt segment inherits the
    transcript's set and reaches parity with the Liftoff path.

    Default ON. Set ``LIFTON_NO_CDS_ATTR_CARRY=1`` to reproduce the pre-change
    ``{Parent}``-only bytes (manuscript reproduction / A-B baseline).
    """
    return os.environ.get("LIFTON_NO_CDS_ATTR_CARRY", "").lower() not in ("1", "true", "yes")


def _descriptive_cds_attrs(attributes):
    """Return a FRESH dict holding only the describing attributes."""
    if not attributes:
        return {}
    return {
        key: value for key, value in attributes.items()
        if key not in _CDS_STRUCTURAL_ATTR_KEYS
    }


def _build_cds_attributes(parent, source_attrs=None, attr_template=None):
    """Build the attribute dict for a rebuilt CDS row.

    Always returns a NEW dict. ``protein_maximization.create_lifton_entries``
    hands every chained CDS of a transcript the SAME ``parent_attrs`` object, so
    reusing it here would let the per-segment ``ID`` assignment in
    ``_normalize_cds_ids`` clobber every sibling. The transcript-level template
    wins over the individual source row so that all segments of one CDS describe
    themselves identically, as the reference does.
    """
    attributes = {}
    if _cds_attr_carry_enabled():
        attributes = _descriptive_cds_attrs(attr_template)
        if not attributes:
            attributes = _descriptive_cds_attrs(source_attrs)
    attributes["Parent"] = parent
    return attributes


class Lifton_ORF:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Lifton_Status:
    def __init__(self):
        self.liftoff = 0
        self.miniprot = 0
        self.lifton_dna = 0
        self.lifton_aa = 0
        self.eval_dna = 0
        self.eval_aa = 0
        self.annotation = None
        self.status = []


class Lifton_Alignment:
    def __init__(self, extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref, cdss_protein_boundary, cdss_protein_aln_boundary, extracted_seq, reference_seq, db_entry):
        self.identity = extracted_identity
        self.cds_children = cds_children
        self.query_aln = alignment_query
        self.comp = alignment_comp
        self.ref_aln = alignment_ref
        self.cdss_protein_boundaries = cdss_protein_boundary
        self.cdss_protein_aln_boundaries = cdss_protein_aln_boundary
        self.query_seq = extracted_seq
        self.ref_seq = reference_seq
        self.db_entry = db_entry

    def write_alignment(self, outdir, tool_name, mutation, trans_id):
        outdir_tool = outdir+"/"+tool_name+"/"+mutation+"/"
        os.makedirs(outdir_tool, exist_ok=True)
        outfile= outdir_tool+trans_id+".fa"
        with open(outfile, "w") as fw:
            fw.write("> Reference\n")
            fw.write(self.ref_aln + "\n")
            fw.write("> Target\n")
            fw.write(self.query_aln + "\n")

    
class Lifton_feature:
    def __init__(self, id):
        self.id = id
        self.copy_num = 0
        self.is_protein_coding = False
        self.is_non_coding = False
        # raw GFF feature type (e.g. gene, mRNA, tRNA, lnc_RNA, pseudogene) and
        # biotype, captured for all-feature-type completeness reporting. Metadata
        # only — never affects the GFF3 written to -o.
        self.feature_type = None
        self.biotype = None
        self.children = set()


class Lifton_GENE:
    def __init__(self, ref_gene_id, gffutil_entry_gene, ref_gene_attrs,
                 tree_dict, ref_features_dict, args, tmp=False,
                 state_journal=None):
        ###########################
        # Assigning the reference gene & attributes
        ###########################
        self.entry = gffutil_entry_gene
        self.entry.source = "LiftOn"
        self.is_protein_coding = False
        self.is_non_coding = False
        self.transcripts = {}
        self.ref_gene_id = ref_gene_id
        if state_journal is None:
            self.copy_num = self.__get_gene_copy(ref_features_dict)
        else:
            # Parallel Step 7 assigns copy IDs through a short ordered gate.
            # The worker records, but does not apply, the corresponding global
            # counter/tree mutations; ordered ``consume`` commits them later.
            self.copy_num = state_journal.allocate_gene_copy(
                ref_gene_id, self.entry, ref_features_dict,
            )
        self.tmp = tmp
        self.entry.attributes = ref_gene_attrs
        # Bug fix #1 (Phase 5): build the ID as a list-of-str per the
        # gffutils attribute contract. Previously a bare string was
        # assigned, and the subsequent `[0]` reduced it to one character.
        gene_id = (f"{self.ref_gene_id}_{self.copy_num}"
                   if self.copy_num > 0 else self.ref_gene_id)
        self.entry.attributes["ID"] = [gene_id]
        if self.copy_num > 0:
            self.entry.attributes["extra_copy_number"] = [str(self.copy_num)]
        if state_journal is None:
            self.__update_gene_copy(ref_features_dict)
        self.entry.id = gene_id
        # V2.8: route via _make_interval so single-base genes don't crash.
        # Only `gene`-type loci join the interval tree that gates miniprot
        # protein-rescue suppression (Step 8). That was always the tree's
        # contract (features was historically ["gene"]); since the default
        # gene-like lift (Iter 12) now also processes pseudogenes/ncRNA_genes,
        # adding them here would let them suppress a miniprot coding-gene rescue
        # just by overlapping it — so they are deliberately kept OUT. On a
        # gene-only ref this is a no-op (byte-identical, 24-cell green).
        if self.entry.featuretype == "gene" and state_journal is None:
            from lifton.intervals import _make_interval
            gene_interval = _make_interval(
                self.entry.start, self.entry.end, self.entry.id,
            )
            if self.entry.seqid not in tree_dict.keys():
                tree_dict[self.entry.seqid] = IntervalTree()
            tree_dict[self.entry.seqid].add(gene_interval)
        # Decide if its type
        # Check for gene type/biotype attribute in order of preference
        gene_type_key = None
        if args.annotation_database.upper() == "REFSEQ":
            # For RefSeq, check gene_biotype first, then biotype as fallback
            if "gene_biotype" in self.entry.attributes.keys():
                gene_type_key = "gene_biotype"
            elif "biotype" in self.entry.attributes.keys():
                gene_type_key = "biotype"
        elif args.annotation_database.upper() == "GENCODE" or args.annotation_database.upper() == "ENSEMBL" or args.annotation_database.upper() == "CHESS":
            # For GENCODE/ENSEMBL/CHESS, check gene_type first, then biotype as fallback
            if "gene_type" in self.entry.attributes.keys():
                gene_type_key = "gene_type"
            elif "biotype" in self.entry.attributes.keys():
                gene_type_key = "biotype"
        else:
            # For other databases, try all possible keys in order
            if "gene_biotype" in self.entry.attributes.keys():
                gene_type_key = "gene_biotype"
            elif "gene_type" in self.entry.attributes.keys():
                gene_type_key = "gene_type"
            elif "biotype" in self.entry.attributes.keys():
                gene_type_key = "biotype"
        
        if gene_type_key is not None:
            if self.entry.attributes[gene_type_key][0] == "protein_coding":
                self.is_protein_coding = True
            elif (self.entry.attributes[gene_type_key][0] == "lncRNA" or self.entry.attributes[gene_type_key][0] == "ncRNA"):
                self.is_non_coding = True
        
    def __get_gene_copy(self, ref_features_dict):
        if 'extra_copy_number' in self.entry.attributes:
            return int(self.entry.attributes['extra_copy_number'][0])
        else:
            if self.ref_gene_id in ref_features_dict.keys():
                return ref_features_dict[self.ref_gene_id].copy_num
            else:
                return 0

    def __update_gene_copy(self, ref_features_dict):
        if self.ref_gene_id in ref_features_dict.keys():
            ref_features_dict[self.ref_gene_id].copy_num += 1

    def update_gene_info(self, chromosome, start, end):
        self.entry.seqid = chromosome
        self.entry.start = start
        self.entry.end = end

    def add_miniprot_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs, ref_features_dict):
        Lifton_trans = Lifton_TRANS(ref_trans_id, self.ref_gene_id, self.entry.id, self.copy_num, gffutil_entry_trans, ref_trans_attrs)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans
    
    def update_trans_info(self, trans_id, chromosome, start, end):
        self.transcripts[trans_id].entry.seqid = chromosome
        self.transcripts[trans_id].entry.start = start
        self.transcripts[trans_id].entry.end = end

    def add_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs):
        Lifton_trans = Lifton_TRANS(ref_trans_id, self.ref_gene_id, self.entry.id, self.copy_num, gffutil_entry_trans, ref_trans_attrs)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans

    def add_feature(self, gffutil_entry_trans):
        Lifton_feature = LiftOn_FEATURE(self.entry.id, gffutil_entry_trans, self.copy_num)
        self.transcripts[Lifton_feature.entry.id] = Lifton_feature
        return Lifton_feature

    def add_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.transcripts[trans_id].add_cds(gffutil_entry_cds)
                            
    def orf_search_protein(self, trans_id, ref_trans_id, fai, ref_proteins, fai_trans, lifton_status, eval_only=False, eval_liftoff_chm13=False):
        if trans_id not in self.transcripts.keys():
            return None, False
        if not eval_liftoff_chm13:
            ref_protein_seq = str(ref_proteins[ref_trans_id]) if ref_trans_id in ref_proteins.keys() else None
            ref_trans_seq = str(fai_trans[ref_trans_id]) if ref_trans_id in fai_trans.keys() else None
        else:
            ref_protein_seq = str(ref_proteins["rna-"+ref_trans_id]) if "rna-"+ref_trans_id in ref_proteins.keys() else None
            ref_trans_seq = str(fai_trans["rna-"+ref_trans_id]) if "rna-"+ref_trans_id in fai_trans.keys() else None
        lifton_aln, good_trans = self.transcripts[trans_id].orf_search_protein(fai, ref_protein_seq, ref_trans_seq, lifton_status, is_non_coding=self.is_non_coding, eval_only=eval_only)
        return lifton_aln, good_trans

    def update_cds_list(self, trans_id, cds_list, optimize=False):
        self.transcripts[trans_id].update_cds_list(cds_list, optimize=optimize)
        self.update_boundaries()

    def add_lifton_gene_status_attrs(self, source):
        self.entry.attributes["source"] = [source]

    def add_lifton_trans_status_attrs(self, trans_id, lifton_status):
        self.transcripts[trans_id].add_lifton_trans_status_attrs(lifton_status)

    def write_entry(
        self,
        fw,
        transcripts_stats_dict,
        cds_id_allocator=None,
    ):
        # Iter-24: normalize parent-child containment (extend exons to cover
        # their CDS; set each transcript + the gene span to the child envelope)
        # right before serialisation, so the emitted GFF3 is always valid.
        # No-op on well-formed genes -> byte-identical default.
        self.normalize_containment(cds_id_allocator=cds_id_allocator)
        # GFF3 serialisation extracted to lifton.io.feature_serializer (Iter 19).
        return feature_serializer.write_gene(self, fw, transcripts_stats_dict)

    def update_boundaries(self):
        for key, trans in self.transcripts.items():
            self.entry.start = trans.entry.start if trans.entry.start < self.entry.start else self.entry.start
            self.entry.end = trans.entry.end if trans.entry.end > self.entry.end else self.entry.end

    def normalize_containment(self, cds_id_allocator=None):
        """Normalize every child's containment, then widen the gene span to
        cover all children. No-op when disabled / well-formed. ``self.transcripts``
        can hold both coding ``Lifton_TRANS`` and gene-like ``LiftOn_FEATURE``
        children (the latter added by ``add_feature`` for pseudogenes / ncRNA
        genes / structured features); both now expose ``normalize_containment``,
        and the ``getattr`` guard keeps a future child type without it from
        aborting the whole-genome write phase."""
        if not _containment_normalize_enabled() or not self.transcripts:
            return
        for trans in self.transcripts.values():
            fn = getattr(trans, "normalize_containment", None)
            if fn is not None:
                if hasattr(trans, "exons"):
                    fn(cds_id_allocator=cds_id_allocator)
                else:
                    fn()
        self.entry.start = min(t.entry.start for t in self.transcripts.values())
        self.entry.end = max(t.entry.end for t in self.transcripts.values())

    def print_gene(self):
        print(self.entry)
        for key, trans in self.transcripts.items():
            trans.print_transcript()


class LiftOn_FEATURE:
    def __init__(self, parent_id, gffutil_entry_feature, copy_num):
        self.entry = gffutil_entry_feature
        self.entry.source = "LiftOn"
        self.copy_num = copy_num
        self.features = {}
        self.entry.attributes["Parent"] = [parent_id]
        self.ref_tran_id = None
        # self.ref_tran_id = self.entry.id
        if int(copy_num) > 0:
            feature_id_base = self.entry.id
            # Liftoff/miniprot may have already appended the extra_copy_number to the ID.
            # E.g. exon-1_1 with extra_copy_number = 1. We must safely strip this before appending our new copy_num.
            if "extra_copy_number" in self.entry.attributes:
                prev_copy = str(self.entry.attributes["extra_copy_number"][0])
                if prev_copy != "0" and prev_copy != "":
                    suffix = f"_{prev_copy}"
                    if feature_id_base.endswith(suffix):
                        feature_id_base = feature_id_base[:-len(suffix)]
            else:
                # Fallback to get_ID_base if extra_copy_number is missing
                feature_id_base = coreutils.get_ID_base(feature_id_base)

            self.entry.id = f"{feature_id_base}_{copy_num}"
            self.entry.attributes["ID"] = [self.entry.id]

    def add_feature(self, gffutil_entry_trans):
        Lifton_feature = LiftOn_FEATURE(self.entry.id, gffutil_entry_trans, self.copy_num)
        self.features[Lifton_feature.entry.id] = Lifton_feature
        return Lifton_feature

    # ── Transcript-bearing delegates ────────────────────────────────────────────
    # A 3-level reference hierarchy (e.g. RefSeq gene -> primary_transcript ->
    # miRNA -> exon) makes `process_liftoff` descend with the INTERMEDIATE feature
    # standing in for the gene. Without these methods that recursion raised
    # AttributeError ('LiftOn_FEATURE' object has no attribute 'add_transcript'),
    # which `process_locus` catches per-locus -- so the transcript, its exons and its
    # CDS were dropped from the output with no diagnostic. Mirror the `Lifton_GENE`
    # API over `self.features`, which `write_entry`/`normalize_containment` already
    # recurse, so a nested transcript is emitted under its real parent.

    def add_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs):
        Lifton_trans = Lifton_TRANS(
            ref_trans_id, getattr(self, "ref_gene_id", None), self.entry.id,
            self.copy_num, gffutil_entry_trans, ref_trans_attrs,
        )
        self.features[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans

    def add_exon(self, trans_id, gffutil_entry_exon):
        self.features[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.features[trans_id].add_cds(gffutil_entry_cds)

    def orf_search_protein(self, trans_id, ref_trans_id, fai, ref_proteins,
                           fai_trans, lifton_status, eval_only=False,
                           eval_liftoff_chm13=False):
        child = self.features.get(trans_id)
        if child is None or not hasattr(child, "orf_search_protein"):
            return None, False
        ref_protein_seq = (str(ref_proteins[ref_trans_id])
                           if ref_trans_id in ref_proteins.keys() else None)
        ref_trans_seq = (str(fai_trans[ref_trans_id])
                         if ref_trans_id in fai_trans.keys() else None)
        return child.orf_search_protein(
            fai, ref_protein_seq, ref_trans_seq, lifton_status,
            is_non_coding=False, eval_only=eval_only,
        )

    def add_lifton_gene_status_attrs(self, source):
        self.entry.attributes["source"] = [source]

    def add_lifton_trans_status_attrs(self, trans_id, lifton_status):
        # Only a Lifton_TRANS child takes the 1-argument form; a nested
        # LiftOn_FEATURE has this same 2-argument signature, so dispatch by type.
        child = self.features.get(trans_id)
        if isinstance(child, Lifton_TRANS):
            child.add_lifton_trans_status_attrs(lifton_status)

    def update_boundaries(self):
        for child in self.features.values():
            if child.entry.start < self.entry.start:
                self.entry.start = child.entry.start
            if child.entry.end > self.entry.end:
                self.entry.end = child.entry.end

    def normalize_containment(self, cds_id_allocator=None):
        """Recursively normalize the generic-feature hierarchy (gene-like
        features: pseudogenes, ncRNA genes, structured mobile elements) so a
        feature always contains its children. Only ever EXTENDS this feature's
        span (never shrinks it), so it cannot drop a child's coverage. No-op when
        disabled or childless. This is the ``LiftOn_FEATURE`` counterpart of the
        ``Lifton_TRANS``/``Lifton_GENE`` containment normalization, so the
        write-funnel call in ``Lifton_GENE.normalize_containment`` does not abort
        the whole-genome write phase on a gene-like-feature child (Iter 24)."""
        if not _containment_normalize_enabled() or not self.features:
            return
        for child in self.features.values():
            fn = getattr(child, "normalize_containment", None)
            if fn is not None:
                fn()
        self.entry.start = min([self.entry.start] + [c.entry.start for c in self.features.values()])
        self.entry.end = max([self.entry.end] + [c.entry.end for c in self.features.values()])

    def write_entry(self, fw):
        # GFF3 serialisation extracted to lifton.io.feature_serializer (Iter 19).
        return feature_serializer.write_feature(self, fw)

    def print_feature(self):
        print(self.entry)
        for key, feature in self.features.items():
            feature.print_feature()


class Lifton_TRANS:
    # Descriptive attributes shared by this transcript's reference CDS rows,
    # captured by `add_cds` and reused whenever the merge or the ORF rescue
    # rebuilds the CDS list. Declared on the class so an instance built without
    # `__init__` (as several tests and the locus proxies do) still resolves it.
    # See `_cds_attr_carry_enabled`.
    _cds_attr_template = None

    def __init__(self, ref_trans_id, ref_gene_id, gene_id, copy_num, gffutil_entry_trans, ref_trans_attrs):
        # Assigning the reference transcripts & attributes
        if int(copy_num) > 0:
            trans_id = f"{ref_trans_id}_{copy_num}" 
        else:
            trans_id = ref_trans_id
        self.entry = gffutil_entry_trans
        self.entry.attributes = ref_trans_attrs
        self.entry.id = trans_id
        self.entry.attributes["ID"] = [trans_id]
        self.entry.source = "LiftOn"
        self.exons = []
        self.exon_dic = {}
        # Get reference ID for the gene & transcript.
        self.ref_gene_id = ref_gene_id
        self.ref_tran_id = ref_trans_id
        if 'Parent' in self.entry.attributes:
            self.entry.attributes['Parent'] = [gene_id]
        if 'transcript_id' in self.entry.attributes:
            self.entry.attributes['transcript_id'] = [self.entry.id]

    def add_lifton_trans_status_attrs(self, lifton_status):
        self.entry.attributes["protein_identity"] = [f"{lifton_status.lifton_aa:.3f}"]
        self.entry.attributes["dna_identity"] = [f"{lifton_status.lifton_dna:.3f}"]
        self.entry.attributes["status"] = [lifton_status.annotation]

    def add_exon(self, gffutil_entry_exon):
        gffutil_entry_exon.attributes['Parent'] = [self.entry.id]
        Lifton_exon = Lifton_EXON(gffutil_entry_exon)
        coreutils.custom_bisect_insert(self.exons, Lifton_exon)

    def add_cds(self, gffutil_entry_cds):
        if self._cds_attr_template is None:
            # First reference CDS wins: every segment of one discontinuous CDS
            # carries the same descriptive attributes, so this is the whole
            # transcript's set.
            template = _descriptive_cds_attrs(
                getattr(gffutil_entry_cds, "attributes", None)
            )
            if template:
                self._cds_attr_template = template
        overlapping = [
            exon for exon in self.exons
            if coreutils.segments_overlap_length(
                (exon.entry.start, exon.entry.end),
                (gffutil_entry_cds.start, gffutil_entry_cds.end))[1]
        ]
        if len(overlapping) > 1:
            # "A CDS lies in exactly one exon" is assumed throughout the model
            # but was never checked, so a CDS spanning an intron was attached to
            # every exon it touched -- emitting the same CDS several times, and
            # (because they shared ONE Feature object) letting a later edit to
            # any copy silently rewrite the others.
            logger.log_warning(
                f"CDS {getattr(gffutil_entry_cds, 'id', '<unnamed>')!r} "
                f"({gffutil_entry_cds.start}-{gffutil_entry_cds.end}) spans "
                f"{len(overlapping)} exons of {self.entry.id}; each will carry "
                "its own copy. The reference model is malformed here."
            )
        for exon in overlapping:
            entry = (gffutil_entry_cds if len(overlapping) == 1
                     else coreutils.clone_feature(gffutil_entry_cds))
            entry.attributes['Parent'] = [self.entry.id]
            exon.add_cds(entry)

    def update_gffutil_entry_trans(self, gffutil_entry_trans):
        for key, atr in gffutil_entry_trans.attributes.items():
            self.entry.attributes[key] = atr

    def mv_exon_idx(self, exon_idx):
        if exon_idx < len(self.exons)-1:
            exon_idx += 1
        return exon_idx

    def update_cds_list(self, cds_list, optimize=False):
        idx_exon_itr = 0
        new_exons = []
        # Reverse CDS list if the strand is "-" => CDSs are in small to large order
        if self.entry.strand == "-":
            cds_list.reverse()
        ########################
        # Guard: empty CDS list
        ########################
        if not cds_list:
            # Remove all CDS from exons but keep exon structure
            for exon in self.exons:
                exon.cds = None
            return
        ########################
        # Case 1: only 1 CDS 
        ########################
        if len(cds_list) == 1:
            only_cds = cds_list[0]
            processed_ovp_exons = False
            ovp_exons = []
            last_exon = None   # <-- track last exon to avoid UnboundLocalError
            while idx_exon_itr < len(self.exons):
                exon = self.exons[idx_exon_itr]
                last_exon = exon
                _, ovp = coreutils.segments_overlap_length(
                    (exon.entry.start, exon.entry.end),
                    (only_cds.entry.start, only_cds.entry.end))

                # |eeeeee| |ccccc|
                if exon.entry.end < only_cds.entry.start:
                    new_exons.append(exon)
                elif ovp:
                    cp_exon = copy.deepcopy(exon)
                    ovp_exons.append(cp_exon)
                # |cccc|  |eeee|
                elif exon.entry.start > only_cds.entry.end:
                    if processed_ovp_exons:
                        # The CDS-bearing exon has ALREADY been emitted; this is a
                        # further downstream exon (3' UTR), so emit it as a plain UTR
                        # exon. This used to be gated on `optimize`, so with
                        # optimize=False (the --legacy-merge path and any direct
                        # caller of Lifton_GENE.update_cds_list, whose default is
                        # False) EVERY downstream exon re-entered the else-branch and
                        # re-emitted another copy of the CDS-bearing exon -- duplicate
                        # exon IDs, a doubled CDS and a wrong protein on any
                        # multi-exon transcript whose CDS spans only the leading
                        # exon(s). The guard is correct for both paths.
                        new_exons.append(exon)
                    else:
                        processed_ovp_exons = True
                        merged_exon = copy.deepcopy(exon)
                        if len(ovp_exons) == 0:
                            # CDS lies entirely before this exon and did not
                            # overlap any earlier exon
                            merged_exon.entry.start = only_cds.entry.start
                            merged_exon.entry.end = only_cds.entry.end
                        else:
                            merged_exon.entry.start = ovp_exons[0].entry.start
                            merged_exon.entry.end = ovp_exons[-1].entry.end
                        merged_exon.add_lifton_cds(only_cds, attr_template=self._cds_attr_template)
                        new_exons.append(merged_exon)
                        new_exons.append(exon)
                idx_exon_itr += 1
            if not processed_ovp_exons:
                # CDS extends to or past the end of all exons
                if last_exon is None:
                    # No exons at all — create a synthetic exon from the CDS
                    import copy as _copy
                    # We cannot create a new entry without a template; skip
                    pass
                else:
                    merged_exon = copy.deepcopy(last_exon)
                    if len(ovp_exons) == 0:
                        merged_exon.entry.start = only_cds.entry.start
                        merged_exon.entry.end = only_cds.entry.end
                    else:
                        merged_exon.entry.start = ovp_exons[0].entry.start
                        merged_exon.entry.end = ovp_exons[-1].entry.end
                    merged_exon.add_lifton_cds(only_cds, attr_template=self._cds_attr_template)
                    new_exons.append(merged_exon)
        ########################
        # Case 2: only 1 exon, and >1 CDSs
        ########################
        elif len(self.exons) == 1:
            first_exon = self.exons[0]
            for cds_idx in range(len(cds_list)):
                exon = copy.deepcopy(first_exon)
                cds = cds_list[cds_idx]
                if cds_idx == 0:
                    if exon.entry.start >= cds.entry.start:
                        exon.entry.start = cds.entry.start
                    exon.entry.end = cds.entry.end
                elif cds_idx == len(cds_list)-1:
                    if exon.entry.end <= cds.entry.end:
                        exon.entry.end = cds.entry.end
                    exon.entry.start = cds.entry.start
                else:
                    exon.entry.end = cds.entry.end
                    exon.entry.start = cds.entry.start
                exon.add_lifton_cds(cds, attr_template=self._cds_attr_template)
                new_exons.append(exon)
        ########################
        # Case 3: multiple CDSs and exons
        ########################
        elif len(cds_list) > 1 and len(self.exons) > 1:
            cds_idx = 0
            exon_idx = 0
            cds = cds_list[cds_idx]
            exon = self.exons[exon_idx]
            ################################################
            # Step 1: Find the first overlapping CDS–exon pair and
            #         emit all exons / CDS that come before it.
            ################################################
            init_head_order = False
            # Determine which starts first: exon or CDS
            # If exon starts before CDS start → exons lead (init_head_order=True)
            # If CDS starts before or at exon start → CDSs lead (False)
            if exon.entry.start < cds.entry.start:
                init_head_order = True
            # (if exon starts >= cds start, init_head_order stays False)

            if init_head_order:
                # Exons lead: scan until we hit an overlap
                while exon_idx < len(self.exons) and cds_idx < len(cds_list):
                    exon_start = self.exons[exon_idx].entry.start
                    exon_end   = self.exons[exon_idx].entry.end
                    cds_start  = cds_list[cds_idx].entry.start
                    cds_end    = cds_list[cds_idx].entry.end
                    _, ovp = coreutils.segments_overlap_length(
                        (exon_start, exon_end), (cds_start, cds_end))
                    # Exon ends before CDS begins → emit exon as UTR
                    if exon_end < cds_start:
                        new_exon = copy.deepcopy(self.exons[exon_idx])
                        new_exon.update_exon_info(exon_start, exon_end)
                        new_exons.append(new_exon)
                        exon_idx += 1
                    # Overlap found → emit first combined exon+CDS, then break
                    elif ovp:
                        new_exon = copy.deepcopy(self.exons[exon_idx])
                        new_exon.update_exon_info(
                            min(exon_start, cds_start), cds_end)
                        new_exon.add_lifton_cds(cds_list[0], attr_template=self._cds_attr_template)
                        new_exons.append(new_exon)
                        exon_idx += 1
                        cds_idx += 1
                        break
                    # CDS ends before exon starts → emit CDS as synthetic exon
                    elif exon_start > cds_end:
                        new_exon = copy.deepcopy(self.exons[exon_idx])
                        new_exon.update_exon_info(cds_start, cds_end)
                        new_exon.add_lifton_cds(cds_list[0], attr_template=self._cds_attr_template)
                        new_exons.append(new_exon)
                        cds_idx += 1
                        break
            else:
                # CDSs lead (or start together): just advance past initial CDSs
                # that precede the first exon boundary, emitting synthetic exons
                while cds_idx < len(cds_list) - 1:
                    cds_start = cds_list[cds_idx].entry.start
                    cds_end   = cds_list[cds_idx].entry.end
                    exon_start = exon.entry.start
                    _, ovp = coreutils.segments_overlap_length(
                        (exon.entry.start, exon.entry.end), (cds_start, cds_end))
                    if ovp:
                        # First overlap found — emit and move on
                        new_exon = copy.deepcopy(exon)
                        new_exon.update_exon_info(cds_start, cds_end)
                        new_exon.add_lifton_cds(cds_list[cds_idx], attr_template=self._cds_attr_template)
                        new_exons.append(new_exon)
                        cds_idx += 1
                        break
                    elif cds_end < exon_start:
                        # CDS comes before first exon
                        new_exon = copy.deepcopy(exon)
                        new_exon.update_exon_info(cds_start, cds_end)
                        new_exon.add_lifton_cds(cds_list[cds_idx], attr_template=self._cds_attr_template)
                        new_exons.append(new_exon)
                        cds_idx += 1
                    else:
                        break
            ################################################
            # Step 2: emit inner CDS blocks (all but the last)
            ################################################
            while cds_idx < len(cds_list) - 1:
                cds = cds_list[cds_idx]
                new_exon = copy.deepcopy(exon)
                new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                new_exon.add_lifton_cds(cds, attr_template=self._cds_attr_template)
                new_exons.append(new_exon)
                cds_idx += 1
            ################################################
            # Step 3: emit the last CDS, merging with its overlapping exon
            ################################################
            last_cds_processed = False
            cds = cds_list[cds_idx]
            cds_start = cds.entry.start
            cds_end   = cds.entry.end
            while exon_idx < len(self.exons):
                exon = self.exons[exon_idx]
                exon_start = exon.entry.start
                exon_end   = exon.entry.end
                _, ovp = coreutils.segments_overlap_length(
                    (exon_start, exon_end), (cds_start, cds_end))
                if exon_end < cds_start:
                    exon_idx += 1
                    continue
                elif ovp:
                    last_cds_processed = True
                    new_exon = copy.deepcopy(exon)
                    new_exon.update_exon_info(cds_start, max(cds_end, exon_end))
                    new_exon.add_lifton_cds(cds, attr_template=self._cds_attr_template)
                    new_exons.append(new_exon)
                    exon_idx += 1
                    break
                elif exon_start > cds_end:
                    last_cds_processed = True
                    new_exon = copy.deepcopy(exon)
                    new_exon.update_exon_info(cds_start, cds_end)
                    new_exon.add_lifton_cds(cds, attr_template=self._cds_attr_template)
                    new_exons.append(new_exon)
                    break
            ################################################
            # Step 4: append any remaining (3' UTR) exons
            ################################################
            while exon_idx < len(self.exons):
                exon = self.exons[exon_idx]
                new_exon = copy.deepcopy(exon)
                new_exon.update_exon_info(new_exon.entry.start, new_exon.entry.end)
                new_exons.append(new_exon)
                exon_idx += 1
            ################################################
            # Step 5: last CDS did not overlap any exon → emit as synthetic exon
            ################################################
            if not last_cds_processed:
                new_exon = copy.deepcopy(self.exons[-1])
                new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                new_exon.add_lifton_cds(cds, attr_template=self._cds_attr_template)
                new_exons.append(new_exon)
        self.exons = new_exons
        if self.exons:  # guard against empty exon list
            if _containment_normalize_enabled():
                # FIX 3 (root): the rebuilt exon list can be non-monotonic — the
                # chaining interleaves miniprot+Liftoff CDS and Case-5 appends a
                # synthetic exon from a deep-copied template — so update_boundaries
                # (exons[0].start..exons[-1].end) could invert the span (start>end)
                # transiently, before the write-funnel normalize_containment fixes
                # it (Iter-24). Sort here so the span is valid AT THE SOURCE for any
                # consumer between here and the writer (the best-of-outcome merge
                # snapshots the span; ORF rescue walks exons). Byte-neutral on the
                # default path (normalize_containment re-sorts to the same order);
                # gated on the same flag so LIFTON_NO_CONTAINMENT_NORMALIZE=1 still
                # reproduces the exact pre-Iteration-24 escape-hatch bytes.
                self.exons.sort(key=lambda e: (e.entry.start, e.entry.end))
            self.update_boundaries()

    def get_coding_seq(self, fai, include_sequence=True):
        """Collect the CDS children, their lengths, and (optionally) the
        spliced coding sequence.

        ``include_sequence=False`` skips the per-CDS ``entry.sequence(fai)``
        fetch and the string concatenation. ``align.LiftOn_translate`` -- the
        only caller, and one that runs up to three times per transcript inside
        the Step-7 hot loop -- immediately OVERWRITES the coding sequence with
        ``get_coding_trans_seq``'s, so every one of those fetches was wasted.
        """
        coding_seq = ""
        cdss_lens = []
        cds_children = []
        for exon in self.exons:
            if exon.cds is not None:
                cds_children.append(coreutils.clone_feature(exon.cds.entry))
                # Chaining the CDS features
                p_seq = exon.cds.entry.sequence(fai) if include_sequence else ""
                if exon.cds.entry.strand == '-':
                    coding_seq = p_seq + coding_seq
                    cdss_lens.insert(0, exon.cds.entry.end - exon.cds.entry.start + 1)
                elif exon.cds.entry.strand == '+':
                    coding_seq = coding_seq + p_seq
                    cdss_lens.append(exon.cds.entry.end - exon.cds.entry.start + 1)
        return coding_seq, cds_children, cdss_lens

    def get_coding_trans_seq(self, fai):
        trans_seq = ""
        coding_seq = ""
        accum_cds_length = 0
        lcl_exons = []
        lcl_exons = self.exons
        if len(self.exons) > 0 and self.exons[0].entry.strand == '-':
            lcl_exons = self.exons[::-1]
        for exon in lcl_exons:
            # Chaining the exon features
            p_trans_seq = exon.entry.sequence(fai)
            p_trans_seq = Seq(p_trans_seq).upper()
            trans_seq = trans_seq + p_trans_seq
            if exon.cds is not None:
                # Updating CDS frame
                exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
                # Chaining the CDS features
                p_seq = exon.cds.entry.sequence(fai)
                coding_seq = coding_seq + p_seq
        if trans_seq != None:
            trans_seq = str(trans_seq).upper()
        if coding_seq != None:
            coding_seq = str(coding_seq).upper()
        return coding_seq, trans_seq

    def translate_coding_seq(self, coding_seq):
        protein_seq = None
        if coding_seq != "":
            # Biopython has historically truncated a terminal partial codon
            # while warning that it may become an error. Preserve that result
            # explicitly so frameshift evaluation remains warning-free.
            complete_length = len(coding_seq) - (len(coding_seq) % 3)
            protein_seq = str(Seq(coding_seq[:complete_length]).translate())
        return protein_seq

    def align_coding_seq(self, protein_seq, ref_protein_seq, lifton_status):
        if ref_protein_seq == "" or ref_protein_seq == None:
            lifton_aa_aln = None
            peps = None
        elif protein_seq == "" or protein_seq == None:
            lifton_aa_aln = None
            peps = None
        else:
            peps = protein_seq.split("*")
            lifton_aa_aln = align.protein_align(protein_seq, ref_protein_seq)
            # Update lifton sequence identity
            lifton_status.lifton_aa = max(lifton_status.lifton_aa, lifton_aa_aln.identity)
        return lifton_aa_aln, peps

    def align_trans_seq(self, trans_seq, ref_trans_seq, lifton_status):
        if ref_trans_seq == "" or ref_trans_seq == None:
            lifton_tran_aln = None
        elif trans_seq == "" or trans_seq == None:
            lifton_tran_aln = None
        else:
            lifton_tran_aln = align.trans_align(trans_seq, ref_trans_seq)
            lifton_status.lifton_dna = lifton_tran_aln.identity
        return lifton_tran_aln

    def orf_search_protein(self, fai, ref_protein_seq, ref_trans_seq, lifton_status, is_non_coding, eval_only=False):
        coding_seq, trans_seq = self.get_coding_trans_seq(fai)
        protein_seq = self.translate_coding_seq(coding_seq)
        # Aligning the LiftOn protein & DNA sequences
        lifton_aa_aln, peps = self.align_coding_seq(protein_seq, ref_protein_seq, lifton_status)
        lifton_tran_aln = self.align_trans_seq(trans_seq, ref_trans_seq, lifton_status)
        # GH #46: the CDS is a contiguous block in the spliced transcript
        # (5'UTR + CDS + 3'UTR), so locate it to scope frameshift detection to the
        # coding region (a UTR indel must not be labelled a coding frameshift).
        cds_span = None
        if coding_seq and trans_seq:
            _cds_start = trans_seq.find(coding_seq)
            if _cds_start >= 0:
                cds_span = (_cds_start, _cds_start + len(coding_seq))
        variants.find_variants(lifton_tran_aln, lifton_aa_aln, lifton_status, peps, is_non_coding, cds_span=cds_span)
        ORF_search = False
        for mutation in lifton_status.status:
            # identical # synonymous 
            # (1) inframe_insertion # (2) inframe_deletion # (3) nonsynonymous 
            # (4) frameshift # (5) start_lost # (6) stop_missing # (7) stop_codon_gain
            # Adding mutations in to entry.attributes
            if mutation != "identical":
                if "mutation" not in self.entry.attributes:
                    self.entry.attributes["mutation"] = [mutation]
                else:
                    self.entry.attributes["mutation"].append(mutation)
            # ORF searching for these four types of mutations
            # frameshift_orf_threshold = 0.8
            # and lifton_aa_aln.identity < frameshift_orf_threshold)
            if mutation == "stop_missing" or mutation == "stop_codon_gain" or mutation == "frameshift"  or mutation == "start_lost":
                ORF_search = True
        if ORF_search and eval_only==False:
            self.__find_orfs(trans_seq, ref_protein_seq, lifton_aa_aln, lifton_status)
        return lifton_tran_aln, lifton_aa_aln

    def __find_orfs(self, trans_seq, ref_protein_seq, lifton_aln, lifton_status):
        """
        Scan the full spliced transcript sequence in all 3 reading frames and
        keep the best ORF per frame (longest ORF that passes the length-growth
        filter), then pick the overall best by protein identity against the
        reference.

        Bug-fixes applied:
          • orf_list now stores exactly one ORF per frame (the longest seen so
            far in that frame), so we don't compare dozens of nested ORFs.
          • The outer ATG scan advances `i` past the end of any ORF it just
            found, avoiding duplicated / overlapping ORF candidates.
          • orf_idx_e is initialised to 0 (not i) before the inner loop so
            the `orf_idx_s < orf_idx_e` guard is correct.
        """
        # Coerce ONCE. The scan below used to wrap every codon slice in str()
        # to tolerate a Bio.Seq; doing it here keeps that safety while letting
        # the inner loop compare plain string slices.
        trans_seq = str(trans_seq).upper()
        start_codon = "ATG"
        stop_codons = {"TAA", "TAG", "TGA"}
        # One best ORF per frame (longest)
        best_orf_per_frame = [None, None, None]  # index = frame (0,1,2)
        max_orf_len = [0, 0, 0]

        for frame in range(3):
            i = frame
            while i < len(trans_seq):
                if trans_seq[i:i + 3] == start_codon:
                    # Scan forward for a stop codon
                    orf_idx_s = i
                    orf_idx_e = 0  # will be set when stop found
                    found_stop = False
                    for j in range(i, len(trans_seq), 3):
                        if trans_seq[j:j + 3] in stop_codons:
                            orf_idx_e = j + 3
                            found_stop = True
                            break
                    if not found_stop:
                        # This scan already covered every codon from `i` to the
                        # end of the sequence IN THIS FRAME, so no later ATG in
                        # the frame can find a stop either. Continuing would
                        # re-scan the same tail from each subsequent ATG, which
                        # is where the worst case turned quadratic.
                        break
                    if orf_idx_s < orf_idx_e:
                        curr_orf_len = orf_idx_e - orf_idx_s
                        if curr_orf_len > max_orf_len[frame]:
                            max_orf_len[frame] = curr_orf_len
                            best_orf_per_frame[frame] = Lifton_ORF(
                                orf_idx_s, orf_idx_e)
                        # Advance past this ORF to avoid nested duplicates
                        i = orf_idx_e
                        continue
                i += 3

        orf_list = [orf for orf in best_orf_per_frame if orf is not None]

        final_orf = None
        max_identity = 0.0
        for orf in orf_list:
            orf_DNA_seq = trans_seq[orf.start:orf.end]
            orf_protein_seq = str(Seq(orf_DNA_seq).translate())
            orf_parasail_res = align.parasail_align_protein_base(
                orf_protein_seq, ref_protein_seq)
            orf_matches, orf_length = get_id_fraction.get_AA_id_fraction(
                orf_parasail_res.traceback.ref,
                orf_parasail_res.traceback.query)
            orf_identity = orf_matches / orf_length if orf_length > 0 else 0.0
            if orf_identity > max_identity:
                max_identity = orf_identity
                final_orf = orf

        # Only update CDS if the new ORF improves identity by at least 1 %
        threshold_orf = 0.01
        if final_orf is not None and max_identity > (lifton_status.lifton_aa + threshold_orf):
            lifton_status.lifton_aa = max_identity
            self.__update_cds_boundary(final_orf)

    def __update_cds_boundary(self, final_orf):
        if self.entry.strand == "+":
            self.__iterate_exons_update_cds(final_orf, self.exons, "+")
        elif self.entry.strand == "-":
            self.__iterate_exons_update_cds(final_orf, self.exons[::-1], "-")            

    def __iterate_exons_update_cds(self, final_orf, exons, strand):
        accum_exon_length = 0
        accum_cds_length = 0
        for exon_idx, exon in enumerate(exons):
            curr_exon_len = exon.entry.end - exon.entry.start + 1
            if accum_exon_length <= final_orf.start:
                if final_orf.start < accum_exon_length+curr_exon_len:
                    # ORF starts in this exon. 
                    # Does it ALSO end in this exon?
                    is_single_exon_orf = (final_orf.end <= accum_exon_length + curr_exon_len)
                    
                    if exon.cds is not None:
                        if strand == "+":
                            exon.cds.entry.start = exon.entry.start + (final_orf.start - accum_exon_length)
                            if is_single_exon_orf:
                                exon.cds.entry.end = exon.entry.start + (final_orf.end - accum_exon_length) - 1
                            else:
                                exon.cds.entry.end = exon.entry.end
                        elif strand == "-":
                            exon.cds.entry.end = exon.entry.end - (final_orf.start - accum_exon_length)
                            if is_single_exon_orf:
                                exon.cds.entry.start = exon.entry.end - (final_orf.end - accum_exon_length) + 1
                            else:
                                exon.cds.entry.start = exon.entry.start 
                    else:
                        if strand == "+":
                            st = exon.entry.start + (final_orf.start - accum_exon_length)
                            en = (exon.entry.start + (final_orf.end - accum_exon_length) - 1) if is_single_exon_orf else exon.entry.end
                            exon.add_novel_lifton_cds(exon.entry, st, en, attr_template=self._cds_attr_template)
                        elif strand == "-":
                            en = exon.entry.end - (final_orf.start - accum_exon_length)
                            st = (exon.entry.end - (final_orf.end - accum_exon_length) + 1) if is_single_exon_orf else exon.entry.start
                            exon.add_novel_lifton_cds(exon.entry, st, en, attr_template=self._cds_attr_template)
                            
                    exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                    accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
                else:
                    # ORF hasn't started yet. No CDS should be created
                    exon.cds = None
            elif final_orf.start < accum_exon_length and accum_exon_length < final_orf.end:
                if final_orf.end <= accum_exon_length+curr_exon_len:
                    # Create the last partial CDS
                    if exon.cds is not None:
                        if strand == "+":
                            exon.cds.entry.start = exon.entry.start
                            exon.cds.entry.end = exon.entry.start + (final_orf.end - accum_exon_length)-1
                        elif strand == "-":
                            exon.cds.entry.end = exon.entry.end
                            exon.cds.entry.start = exon.entry.end - (final_orf.end - accum_exon_length) + 1
                    else:
                        if strand == "+":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.start + (final_orf.end - accum_exon_length)-1, attr_template=self._cds_attr_template)
                        elif strand == "-":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.end - (final_orf.end - accum_exon_length) + 1, exon.entry.end, attr_template=self._cds_attr_template)
                else:
                    # Keep the original full CDS / extend the CDS to full exon length
                    if exon.cds is None:
                        exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.end, attr_template=self._cds_attr_template)
                    else:
                        exon.cds.update_CDS_info(exon.entry.start, exon.entry.end)
                exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
            elif final_orf.end <= accum_exon_length:
                # ORF has already ended. No CDS should be created
                exon.cds = None
            accum_exon_length += curr_exon_len

    def __get_cds_frame(self, accum_cds_length):
        return (3 - accum_cds_length%3)%3

    def write_entry(self, fw):
        # GFF3 serialisation extracted to lifton.io.feature_serializer (Iter 19).
        return feature_serializer.write_trans(self, fw)

    def update_boundaries(self):
        self.entry.start = self.exons[0].entry.start
        self.entry.end = self.exons[-1].entry.end

    def normalize_containment(self, cds_id_allocator=None):
        """Guarantee GFF3 parent-child containment for this transcript.

        ORF-boundary patching / chaining can leave a CDS extending a few bp
        past its exon (e.g. capturing a stop codon on a frameshift-corrected
        RefSeq model) or rebuild the exons out of coordinate order, yielding
        child features outside the mRNA span (invalid GFF3). A CDS must lie
        within a transcribed exon, so extend the exon to cover its CDS; then
        sort the exons and set the transcript span to the exon envelope. A
        no-op on well-formed transcripts (exons already contain their CDS, are
        already sorted, and the span already equals the envelope), so the
        default output is byte-identical except on the previously-invalid
        transcripts. Disable with LIFTON_NO_CONTAINMENT_NORMALIZE=1.
        """
        if not _containment_normalize_enabled() or not self.exons:
            return
        for exon in self.exons:
            if exon.cds is not None:
                if exon.cds.entry.start < exon.entry.start:
                    exon.entry.start = exon.cds.entry.start
                if exon.cds.entry.end > exon.entry.end:
                    exon.entry.end = exon.cds.entry.end
        self.exons.sort(key=lambda e: (e.entry.start, e.entry.end))
        self.entry.start = min(e.entry.start for e in self.exons)
        self.entry.end = max(e.entry.end for e in self.exons)
        # Ensure exon IDs are unique within the transcript. The chaining can
        # rebuild the exon structure while reusing stale reference "exon-<acc>-N"
        # IDs, producing DUPLICATE exon IDs (invalid GFF3 — the dominant validator
        # error). Re-number 5'->3' from the shared base, but ONLY when a real
        # collision exists, so well-formed transcripts (already-unique, strand-
        # ordered IDs) are left untouched / byte-identical.
        exon_ids = [e.entry.id for e in self.exons]
        if len(set(exon_ids)) != len(exon_ids):
            base = None
            for e in self.exons:
                m = re.match(r"^(.*)-\d+$", e.entry.id or "")
                if m:
                    base = m.group(1)
                    break
            if base is None:
                base = "exon-" + (self.entry.id or "")
            ordered = sorted(self.exons, key=lambda e: e.entry.start,
                             reverse=(self.entry.strand == "-"))
            for i, exon in enumerate(ordered, 1):
                new_id = f"{base}-{i}"
                exon.entry.id = new_id
                exon.entry.attributes["ID"] = [new_id]
        # GH #32/#8: give every CDS an ID. Preserve an existing stable source ID
        # (for example NCBI's cds-NP_*) and use it to fill missing segments.
        # Chaining/miniprot can rebuild every segment without an ID; only then
        # synthesize the shared discontinuous-feature ID from the transcript.
        cds_features = [
            exon.cds for exon in self.exons if exon.cds is not None
        ]
        if not cds_features:
            return

        def first_id(values):
            if isinstance(values, str):
                values = [values]
            return next(
                (str(value) for value in values if value not in (None, "")),
                None,
            )

        reserved_ids = {
            feature_id
            for feature_id in (
                first_id(exon.entry.attributes.get("ID", []))
                for exon in self.exons
            )
            if feature_id
        }
        if self.entry.id:
            reserved_ids.add(self.entry.id)

        ref_trans_id = getattr(self, "ref_tran_id", None)
        copy_suffix = ""
        if (
            ref_trans_id
            and self.entry.id
            and self.entry.id.startswith(f"{ref_trans_id}_")
        ):
            copy_suffix = self.entry.id[len(ref_trans_id):]
            if not copy_suffix[1:].isdigit():
                copy_suffix = ""

        def preservable_cds_id(cds):
            cds_id = getattr(cds, "_source_id_base", None)
            if not cds_id or cds_id in reserved_ids:
                return None
            if copy_suffix:
                cds_id += copy_suffix
            return None if cds_id in reserved_ids else cds_id

        normalized_ids = [
            preservable_cds_id(cds) for cds in cds_features
        ]
        allocator = cds_id_allocator or CdsIdAllocator()
        parent_id = str(self.entry.id or "")
        source_ids = [
            getattr(cds, "_source_id_base", None)
            for cds in cds_features
        ]
        for cds_id in dict.fromkeys(
            identifier for identifier in normalized_ids if identifier
        ):
            source_id = next(
                source
                for identifier, source in zip(normalized_ids, source_ids)
                if identifier == cds_id
            )
            if allocator.reserve_explicit(
                cds_id,
                parent_id,
                source_identifier=source_id,
            ):
                continue
            normalized_ids = [
                None if identifier == cds_id else identifier
                for identifier in normalized_ids
            ]
        existing_ids = list(dict.fromkeys(
            cds_id for cds_id in normalized_ids if cds_id
        ))
        synthetic_base = (
            "cds-" + re.sub(r"^rna-", "", (self.entry.id or ""))
        )
        needs_synthetic = (
            not existing_ids
            or (len(existing_ids) > 1 and any(
                identifier is None for identifier in normalized_ids
            ))
        )
        synthetic_id = (
            allocator.allocate_synthetic(
                synthetic_base,
                parent_id,
                local_reserved_ids=reserved_ids | set(existing_ids),
            )
            if needs_synthetic else None
        )

        if len(existing_ids) <= 1:
            shared_id = existing_ids[0] if existing_ids else synthetic_id
            for cds in cds_features:
                cds.entry.attributes["ID"] = [shared_id]
            return

        # Multiple pre-existing IDs can be valid segment-level identifiers.
        # Do not collapse them; fill only gaps. Copies still need distinct IDs
        # from the source transcript, so append the transcript's copy suffix.
        for cds, cds_id in zip(cds_features, normalized_ids):
            cds.entry.attributes["ID"] = [cds_id or synthetic_id]

    def print_transcript(self):
        print(f"\t{self.entry}")
        for exon in self.exons:
            exon.print_exon()


class Lifton_EXON:
    def __init__(self, gffutil_entry_exon):
        gffutil_entry_exon.source = "LiftOn"
        gffutil_entry_exon.featuretype = "exon"
        self.entry = gffutil_entry_exon
        if 'extra_copy_number' in self.entry.attributes: self.entry.attributes.pop('extra_copy_number')
        self.cds = None

    def __deepcopy__(self, memo):
        # `copy.deepcopy(trans.exons)` in `run_liftoff._snapshot_merge_state` and
        # `copy.deepcopy(exon)` in `update_cds_list` were the two hottest copy
        # sites in the Step-7 profile. An exon holds exactly one gffutils
        # Feature and at most one Lifton_CDS, so a purpose-built clone is both
        # complete and ~6.6x cheaper per Feature. `memo` is honoured so a shared
        # reference stays shared, exactly as deepcopy guarantees.
        existing = memo.get(id(self))
        if existing is not None:
            return existing
        clone = Lifton_EXON.__new__(Lifton_EXON)
        memo[id(self)] = clone
        clone.entry = coreutils.clone_feature(self.entry)
        clone.cds = copy.deepcopy(self.cds, memo) if self.cds is not None else None
        return clone

    def update_exon_info(self, start, end):
        self.cds = None
        self.entry.source = "LiftOn"
        self.entry.start = start
        self.entry.end = end

    def add_cds(self, gffutil_entry_cds):
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        self.cds = Lifton_cds

    def add_novel_lifton_cds(self, gffutil_entry_exon, start, end,
                             attr_template=None):
        gffutil_entry_cds = coreutils.clone_feature(gffutil_entry_exon)
        gffutil_entry_cds.featuretype = "CDS"
        gffutil_entry_cds.start = start
        gffutil_entry_cds.end = end
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        Lifton_cds._source_id = None
        Lifton_cds._source_id_base = None
        # The template, not the exon's own attributes: an exon row carries
        # `gbkey=mRNA` and no `protein_id`, so copying it would mislabel the CDS.
        Lifton_cds.entry.attributes = _build_cds_attributes(
            self.entry.attributes['Parent'], attr_template=attr_template,
        )
        self.cds = Lifton_cds

    def add_lifton_cds(self, Lifton_cds, attr_template=None):
        if Lifton_cds is not None:
            # Source ID/copy provenance lives on ``Lifton_CDS`` until the
            # enabled write-funnel normalizer validates and materializes it, so
            # ``ID`` is never carried here. The describing attributes ARE, so a
            # rebuilt CDS reads like the reference row it came from.
            Lifton_cds.entry.attributes = _build_cds_attributes(
                self.entry.attributes['Parent'],
                source_attrs=Lifton_cds.entry.attributes,
                attr_template=attr_template,
            )
        self.cds = Lifton_cds

    def write_entry(self, fw):
        # GFF3 serialisation extracted to lifton.io.feature_serializer (Iter 19).
        return feature_serializer.write_exon(self, fw)

    def print_exon(self):
        print(f"\t\t{self.entry}")
        if self.cds != None:
            self.cds.print_cds()


class Lifton_CDS:
    def __init__(self, gffutil_entry_cds):
        gffutil_entry_cds.source = "LiftOn"
        gffutil_entry_cds.featuretype = "CDS"
        self.entry = gffutil_entry_cds
        source_id = self.entry.attributes.get("ID", [])
        if isinstance(source_id, str):
            self._source_id = source_id or None
        else:
            self._source_id = next(
                (
                    str(identifier)
                    for identifier in source_id
                    if identifier not in (None, "")
                ),
                None,
            )
        source_copy_number = self.entry.attributes.get(
            "extra_copy_number", []
        )
        if isinstance(source_copy_number, str):
            source_copy_number = [source_copy_number]
        self._source_copy_number = (
            str(source_copy_number[0]) if source_copy_number else None
        )
        copy_suffix = (
            f"_{self._source_copy_number}"
            if self._source_copy_number not in (None, "", "0") else ""
        )
        self._source_id_base = self._source_id
        if (
            self._source_id_base
            and copy_suffix
            and self._source_id_base.endswith(copy_suffix)
        ):
            self._source_id_base = self._source_id_base[:-len(copy_suffix)]
        if 'extra_copy_number' in self.entry.attributes: self.entry.attributes.pop('extra_copy_number')

    def __deepcopy__(self, memo):
        # Copied verbatim rather than re-run through __init__: the source-ID
        # provenance (_source_id / _source_id_base / _source_copy_number) is
        # derived from attributes that __init__ POPS, so re-deriving it from a
        # clone would silently lose the copy suffix. See Lifton_EXON.__deepcopy__.
        existing = memo.get(id(self))
        if existing is not None:
            return existing
        clone = Lifton_CDS.__new__(Lifton_CDS)
        memo[id(self)] = clone
        clone.__dict__.update(self.__dict__)
        clone.entry = coreutils.clone_feature(self.entry)
        return clone

    def update_CDS_info(self, start, end):
        self.entry.source = "LiftOn"
        self.entry.start = start
        self.entry.end = end

    def write_entry(self, fw):
        # GFF3 serialisation extracted to lifton.io.feature_serializer (Iter 19).
        return feature_serializer.write_cds(self, fw)

    def print_cds(self):
        print(f"\t\t{self.entry}")
