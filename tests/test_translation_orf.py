"""Phase 4.5 Step 2 — translation, ORF rescue, and variant classification.

Drives Lifton_TRANS.translate_coding_seq, __find_orfs,
__update_cds_boundary, and lifton.variants.find_variants through
synthetic ORFs covering each NCBI-relevant edge case.
"""

from __future__ import annotations

import copy

import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from lifton import (
    annotation,
    extract_sequence,
    lifton_class,
    lifton_utils,
    variants,
)


# ---------------------------------------------------------------------------
# Pure translation behaviour (BioPython contract)
# ---------------------------------------------------------------------------

class TestTranslateCodingSeq:
    def _trans(self):
        # Bare Lifton_TRANS-like object — only translate_coding_seq is
        # called, which is a static-ish method on the instance.
        from types import SimpleNamespace
        instance = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
        return instance

    def test_canonical_protein(self):
        coding = "ATG" + "GCT" * 32 + "TAA"
        protein = self._trans().translate_coding_seq(coding)
        assert protein.startswith("M")
        assert protein.endswith("*")
        assert protein.count("*") == 1

    def test_premature_stop_truncates_protein(self):
        # Canonical start, internal stop after 3 codons, then continues
        coding = "ATG" + "GCT" + "TAA" + "GCT" * 30 + "TAA"
        protein = self._trans().translate_coding_seq(coding)
        # BioPython's translate emits both stops; LiftOn does NOT truncate
        # at the first one. Pin current behaviour.
        assert protein.count("*") == 2
        assert protein.startswith("MA*")

    def test_missing_stop_returns_protein_without_terminator(self):
        coding = "ATG" + "GCT" * 65   # 198 bp, no stop codon at all
        protein = self._trans().translate_coding_seq(coding)
        assert protein.startswith("M")
        assert "*" not in protein

    def test_noncanonical_start_does_not_yield_M(self):
        # CTG is an alternative start in some organisms but standard table
        # translates it to L (leucine). LiftOn uses the standard table.
        coding = "CTG" + "GCT" * 32 + "TAA"
        protein = self._trans().translate_coding_seq(coding)
        assert protein.startswith("L")
        assert not lifton_utils.check_protein_valid(protein)

    def test_empty_coding_returns_none(self):
        assert self._trans().translate_coding_seq("") is None

    def test_selenocysteine_tga_treated_as_stop(self):
        """NCBI uses transl_except to mark TGA → Sec/U. LiftOn's plain
        BioPython translate sees TGA as `*`. This pins the legacy gap;
        Phase 4 Step 2 validator should warn when the input GFF carries
        a transl_except attribute that LiftOn cannot honour."""
        coding = "ATG" + "GCT" * 5 + "TGA" + "GCT" * 5 + "TAA"
        protein = self._trans().translate_coding_seq(coding)
        assert protein.count("*") == 2  # internal TGA + terminal TAA


# ---------------------------------------------------------------------------
# Reverse-strand path through extract_sequence
# ---------------------------------------------------------------------------

class TestReverseStrandExtraction:
    def test_reverse_strand_yields_canonical_protein(self,
                                                     fasta_reverse_strand):
        from types import SimpleNamespace
        fa = Fasta(str(fasta_reverse_strand))
        # On the - strand, the ORF lives at chr1:101..202 in fwd coords;
        # extract_sequence reverse-complements when strand is "-".
        # Build a class-like stub for child intervals.
        class _Iv:
            def __init__(self, start, end):
                self.start = start; self.end = end

        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="-"), fa,
            [_Iv(101, 202)],
        )
        protein = str(Seq(seq).translate())
        assert protein.startswith("M")
        assert protein.endswith("*")


# ---------------------------------------------------------------------------
# get_padding_length on frameshift remnants
# ---------------------------------------------------------------------------

class TestPaddingFrameshift:
    @pytest.mark.parametrize("L,pad", [
        (100, 2),  # 100 = 3*33 + 1 -> needs 2 to round up to 102
        (101, 1),  # 101 = 3*33 + 2 -> needs 1
        (99, 0),   # 99 already multiple of 3
    ])
    def test_padding_for_frameshifted_lengths(self, L, pad):
        assert extract_sequence.get_padding_length(L) == pad


# ---------------------------------------------------------------------------
# Variant classification (lifton.variants)
# ---------------------------------------------------------------------------

class TestFindVariants:
    def _aln(self, query, ref, identity=1.0, query_seq=None):
        return lifton_class.Lifton_Alignment(
            extracted_identity=identity,
            cds_children=None,
            alignment_query=query,
            alignment_comp=None,
            alignment_ref=ref,
            cdss_protein_boundary=None,
            cdss_protein_aln_boundary=None,
            extracted_seq=query_seq if query_seq is not None else query,
            reference_seq=ref,
            db_entry=None,
        )

    def test_identical_alignment(self):
        status = lifton_class.Lifton_Status()
        dna = self._aln("ATGGCT", "ATGGCT", identity=1.0)
        prot = self._aln("MA*", "MA*", identity=1.0)
        variants.find_variants(dna, prot, status, ["MA", ""], False)
        assert status.status == ["identical"]

    def test_synonymous_when_protein_identical_dna_diverges(self):
        status = lifton_class.Lifton_Status()
        dna = self._aln("ATGGCC", "ATGGCT", identity=0.83)
        prot = self._aln("MA*", "MA*", identity=1.0)
        variants.find_variants(dna, prot, status, ["MA", ""], False)
        assert status.status == ["synonymous"]

    def test_frameshift_detected_via_query_alignment(self):
        status = lifton_class.Lifton_Status()
        # 1-bp insertion in query (gap of length 1 in ref -> frameshift)
        dna = self._aln("ATG-GCTTAA", "ATGAGCTTAA", identity=0.8)
        prot = self._aln("M*", "MA*", identity=0.5,
                         query_seq="M*")
        variants.find_variants(dna, prot, status, ["M", ""], False)
        assert "frameshift" in status.status

    def test_inframe_insertion_when_3bp_gap_in_ref(self):
        status = lifton_class.Lifton_Status()
        # 3-bp gap in ref -> inframe insertion (not frameshift)
        dna = self._aln("ATGGCTTAA", "ATG---TAA", identity=0.66)
        prot = self._aln("MA*", "M-*", identity=0.66,
                         query_seq="MA*")
        variants.find_variants(dna, prot, status, ["MA", ""], False)
        assert "inframe_insertion" in status.status
        assert "frameshift" not in status.status

    def test_stop_missing_when_peps_has_one_element(self):
        status = lifton_class.Lifton_Status()
        dna = self._aln("ATGGCT", "ATGGCT", identity=0.9)
        prot = self._aln("MA", "MA*", identity=0.5, query_seq="MA")
        variants.find_variants(dna, prot, status, ["MA"], False)
        assert "stop_missing" in status.status

    def test_stop_codon_gain_when_multiple_peps(self):
        status = lifton_class.Lifton_Status()
        dna = self._aln("ATGTAAGCT", "ATGGCTGCT", identity=0.5)
        prot = self._aln("M*A", "MAA", identity=0.33,
                         query_seq="M*A")
        variants.find_variants(dna, prot, status, ["M", "A"], False)
        assert "stop_codon_gain" in status.status

    def test_non_coding_short_circuits(self):
        status = lifton_class.Lifton_Status()
        variants.find_variants(None, None, status, [], True)
        assert status.status == ["non_coding"]

    def test_full_transcript_loss_when_dna_align_none(self):
        status = lifton_class.Lifton_Status()
        variants.find_variants(None, "anything", status, [], False)
        assert status.status == ["full_transcript_loss"]

    def test_no_protein_when_protein_align_none(self):
        status = lifton_class.Lifton_Status()
        dna = self._aln("ATG", "ATG", identity=1.0)
        # identity=1.0 short-circuits to 'identical', so use 0.5
        dna = self._aln("ATG", "ATC", identity=0.66)
        variants.find_variants(dna, None, status, [], False)
        assert status.status == ["no_protein"]


# ---------------------------------------------------------------------------
# is_frameshift / has_stop_codon helpers
# ---------------------------------------------------------------------------

class TestVariantHelpers:
    def test_is_frameshift_true_for_one_bp_gap(self):
        assert variants.is_frameshift("ATG-GCT") is True

    def test_is_frameshift_true_for_two_bp_gap(self):
        assert variants.is_frameshift("ATG--GCT") is True

    def test_is_frameshift_false_for_three_bp_gap(self):
        assert variants.is_frameshift("ATG---GCT") is False

    def test_is_frameshift_false_for_no_gaps(self):
        assert variants.is_frameshift("ATGGCT") is False

    def test_is_frameshift_handles_trailing_gap(self):
        assert variants.is_frameshift("ATGGCT-") is True
        assert variants.is_frameshift("ATGGCT---") is False

    def test_has_stop_codon_true_when_target_has_extra_stop(self):
        assert variants.has_stop_codon("MAGT*", "MA*GT*") is True

    def test_has_stop_codon_false_when_aligned_stops_match(self):
        assert variants.has_stop_codon("MAGT*", "MAGT*") is False


# ---------------------------------------------------------------------------
# ORF rescue end-to-end via Lifton_TRANS.orf_search_protein
# ---------------------------------------------------------------------------

def _build_trans_with_premature_stop(gff_standard, fasta_with_premature_stop,
                                     fake_args, ref_features_dict_one_gene):
    db = annotation.Annotation(
        str(gff_standard), False, False, "create_unique", None, True, False,
    ).db_connection
    gene_entry = copy.deepcopy(db["gene1"])
    gene = lifton_class.Lifton_GENE(
        ref_gene_id="gene1",
        gffutil_entry_gene=gene_entry,
        ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["protein_coding"]},
        tree_dict={},
        ref_features_dict=ref_features_dict_one_gene,
        args=fake_args,
    )
    trans = gene.add_transcript(
        "tx1", copy.deepcopy(db["tx1"]),
        {"ID": ["tx1"], "Parent": ["gene1"]},
    )
    trans.add_exon(copy.deepcopy(db["exon1"]))
    trans.add_exon(copy.deepcopy(db["exon2"]))
    trans.add_cds(copy.deepcopy(db["cds1"]))
    trans.add_cds(copy.deepcopy(db["cds2"]))
    return gene, trans, Fasta(str(fasta_with_premature_stop))


class TestOrfRescueEndToEnd:
    def test_premature_stop_triggers_mutation_attribute(
            self, gff_standard, fasta_with_premature_stop, fake_args,
            ref_features_dict_one_gene):
        gene, trans, fa = _build_trans_with_premature_stop(
            gff_standard, fasta_with_premature_stop, fake_args,
            ref_features_dict_one_gene,
        )
        # Reference protein is the canonical fasta_standard ORF
        from tests.conftest import _build_chrom_with_gene  # noqa
        canon_seq = "ATG" + "GCT" * 32 + "GCT" * 32 + "TAA"
        ref_protein = str(Seq(canon_seq).translate())
        status = lifton_class.Lifton_Status()
        trans.orf_search_protein(
            fa, ref_protein, canon_seq, status,
            is_non_coding=False, eval_only=False,
        )
        # The premature-stop transcript must be flagged
        assert status.status != ["identical"]
        # And `mutation` attribute populated on entry
        assert "mutation" in trans.entry.attributes

    def test_eval_only_skips_orf_rescue(
            self, gff_standard, fasta_with_premature_stop, fake_args,
            ref_features_dict_one_gene):
        gene, trans, fa = _build_trans_with_premature_stop(
            gff_standard, fasta_with_premature_stop, fake_args,
            ref_features_dict_one_gene,
        )
        canon_seq = "ATG" + "GCT" * 64 + "TAA"
        ref_protein = str(Seq(canon_seq).translate())
        status = lifton_class.Lifton_Status()
        # Capture exon coordinates before
        exons_before = [(e.entry.start, e.entry.end) for e in trans.exons]
        trans.orf_search_protein(
            fa, ref_protein, canon_seq, status,
            is_non_coding=False, eval_only=True,   # <-- must NOT mutate CDS
        )
        exons_after = [(e.entry.start, e.entry.end) for e in trans.exons]
        # eval_only=True ⇒ ORF rescue skipped, exon boundaries untouched
        assert exons_before == exons_after
