"""Shared pytest fixtures for the LiftOn legacy-behaviour test suite.

Fixtures are intentionally small and hermetic: each FASTA chromosome is
~600 bp, hand-crafted so the embedded gene yields a syntactically valid
ATG..TAA ORF on the forward strand. No external binaries are invoked
from any test that uses these fixtures.
"""

from __future__ import annotations

import os
import shutil
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# ---------------------------------------------------------------------------
# Sequence builders
# ---------------------------------------------------------------------------

def _wrap(seq: str, width: int = 60) -> str:
    return "\n".join(textwrap.wrap(seq, width)) + "\n"


def _padded(seq: str, total: int) -> str:
    """Right-pad a seed sequence with deterministic 'A' bases up to length."""
    if len(seq) >= total:
        return seq[:total]
    return seq + "A" * (total - len(seq))


def _build_chrom_with_gene(chrom_len: int = 600) -> str:
    """Construct a deterministic chromosome whose 1-based coordinates
    101..199 contain exon1 (ATG..) and 301..399 contain exon2 (..TAA).

    The two exons concatenated form an in-frame ORF whose translation is
    the reference protein used in golden-path expectations.
    """
    base = ["A"] * chrom_len
    # exon1: positions 101..199 (1-based inclusive) -> indices 100..198
    exon1 = "ATG" + "GCT" * 32  # 3 + 96 = 99 nt, no stop codons
    for i, ch in enumerate(exon1):
        base[100 + i] = ch
    # intron 200..300 stays "A"
    # exon2: positions 301..399 -> indices 300..398
    exon2 = ("GCT" * 32) + "TAA"  # 96 + 3 = 99 nt, ends with stop
    for i, ch in enumerate(exon2):
        base[300 + i] = ch
    return "".join(base)


# ---------------------------------------------------------------------------
# File-on-disk fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fasta_standard(tmp_path: Path) -> Path:
    """Single-chromosome FASTA covering the standard gene used below."""
    fp = tmp_path / "ref.fa"
    fp.write_text(">chr1\n" + _wrap(_build_chrom_with_gene()))
    return fp


@pytest.fixture
def fasta_two_chrom(tmp_path: Path) -> Path:
    fp = tmp_path / "ref_two.fa"
    fp.write_text(
        ">chr1\n" + _wrap(_build_chrom_with_gene())
        + ">chr2\n" + _wrap(_padded("ACGTACGT", 600))
    )
    return fp


@pytest.fixture
def fasta_missing_chrom(tmp_path: Path) -> Path:
    """FASTA that does NOT contain chr1; used to exercise the missing-chrom
    early return in extract_sequence.get_dna_sequence."""
    fp = tmp_path / "ref_missing.fa"
    fp.write_text(">chrZ\n" + _wrap(_padded("ACGT", 200)))
    return fp


@pytest.fixture
def gff_standard(tmp_path: Path) -> Path:
    """A minimal but well-formed GFF3 with one protein-coding gene,
    one mRNA, two exons (101..199, 301..399) and matching CDSs."""
    fp = tmp_path / "ref.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\ttest\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\ttest\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )
    return fp


@pytest.fixture
def gff_single_cds(tmp_path: Path) -> Path:
    """Single-exon, single-CDS gene used by Lifton_TRANS update_cds_list
    Case 1 unit tests."""
    fp = tmp_path / "ref_single.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=g_s;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t199\t.\t+\t.\tID=tx_s;Parent=g_s\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=ex_s;Parent=tx_s\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds_s;Parent=tx_s\n"
    )
    return fp


@pytest.fixture
def gff_noncoding(tmp_path: Path) -> Path:
    """ncRNA gene (no CDS) for the noncoding partition path."""
    fp = tmp_path / "ref_nc.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=ncg;gene_biotype=lncRNA\n"
        "chr1\ttest\tlnc_RNA\t101\t199\t.\t+\t.\tID=nctx;Parent=ncg\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=ncex;Parent=nctx\n"
    )
    return fp


@pytest.fixture
def gff_malformed_coords(tmp_path: Path) -> Path:
    """GFF3 whose end < start on one feature. gffutils tolerates this on
    parse but downstream coordinate math should not crash on it."""
    fp = tmp_path / "ref_bad.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t200\t100\t.\t+\t.\tID=gbad;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t200\t100\t.\t+\t.\tID=txbad;Parent=gbad\n"
    )
    return fp


# ---------------------------------------------------------------------------
# In-memory helpers (no I/O)
# ---------------------------------------------------------------------------

@pytest.fixture
def make_gffutils_feature():
    """Factory returning a gffutils.Feature from a single GFF3 line.
    Used by Lifton_TRANS / Lifton_EXON / Lifton_CDS unit tests so they
    can construct entries without spinning up a full FeatureDB."""
    import gffutils

    def _make(seqid="chr1", source="test", featuretype="exon",
              start=1, end=100, strand="+", frame=".", attributes=None,
              feature_id=None):
        attrs = dict(attributes or {})
        if feature_id is not None and "ID" not in attrs:
            attrs["ID"] = [feature_id]
        # Derive feature.id from ID attribute if present (mirrors gffutils
        # parser behaviour). gffutils.Feature does not do this automatically
        # when constructed in-memory.
        derived_id = feature_id
        if derived_id is None and "ID" in attrs and attrs["ID"]:
            derived_id = attrs["ID"][0]
        return gffutils.Feature(
            seqid=seqid, source=source, featuretype=featuretype,
            start=start, end=end, strand=strand, frame=frame,
            attributes=attrs, id=derived_id or "",
        )

    return _make


@pytest.fixture
def fake_args():
    """A SimpleNamespace populated with the args fields that the
    constructors of Lifton_GENE/TRANS read. Default annotation_database
    is RefSeq so gene_biotype is the discriminator."""
    return SimpleNamespace(
        annotation_database="RefSeq",
        evaluation=False,
        evaluation_liftoff_chm13=False,
        debug=False,
    )


@pytest.fixture
def empty_tree_dict():
    return {}


@pytest.fixture
def ref_features_dict_one_gene():
    """Mirrors what get_ref_liffover_features would produce for gene1."""
    from lifton import lifton_class

    feature = lifton_class.Lifton_feature("gene1")
    feature.is_protein_coding = True
    return {"gene1": feature}


# ---------------------------------------------------------------------------
# Cleanup of gffutils SQLite caches across tests
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Phase 4.5 — additional fixtures for boundary / ORF / multi-isoform / dirty
# input tests. Keep above the autouse cleanup fixture.
# ---------------------------------------------------------------------------

@pytest.fixture
def fasta_short_chrom(tmp_path: Path) -> Path:
    """Tiny 50-bp chromosome for boundary-edge tests (1-bp and end-of-chrom
    features). Sequence is `ATG` followed by 47 bases of `A` so the first
    three bases are a start codon."""
    fp = tmp_path / "tiny.fa"
    body = "ATG" + "A" * 47
    fp.write_text(">chr1\n" + _wrap(body))
    return fp


@pytest.fixture
def fasta_with_premature_stop(tmp_path: Path) -> Path:
    """600 bp chromosome whose embedded gene CDS contains an in-frame
    premature TAA at the start of exon 2. exon1 ends after ATG + 96 bp,
    exon2 starts with TAA so the merged CDS translation truncates early."""
    base = ["A"] * 600
    exon1 = "ATG" + "GCT" * 32                     # 99 bp, no stop
    exon2 = "TAA" + ("GCT" * 31) + "GCG" + "TAA"   # 99 bp, internal stop at base 0
    for i, ch in enumerate(exon1):
        base[100 + i] = ch
    for i, ch in enumerate(exon2):
        base[300 + i] = ch
    fp = tmp_path / "premature_stop.fa"
    fp.write_text(">chr1\n" + _wrap("".join(base)))
    return fp


@pytest.fixture
def fasta_no_stop(tmp_path: Path) -> Path:
    """600 bp chromosome whose CDS lacks any of TAA/TAG/TGA."""
    base = ["A"] * 600
    cds = "ATG" + "GCT" * 65   # 198 bp = 66 codons of M+A's, no stop
    for i, ch in enumerate(cds):
        base[100 + i] = ch
    fp = tmp_path / "no_stop.fa"
    fp.write_text(">chr1\n" + _wrap("".join(base)))
    return fp


@pytest.fixture
def fasta_noncanonical_start(tmp_path: Path) -> Path:
    """600 bp chromosome whose CDS starts with CTG (alternative start)."""
    base = ["A"] * 600
    cds = "CTG" + "GCT" * 32 + "TAA"   # 102 bp
    for i, ch in enumerate(cds):
        base[100 + i] = ch
    fp = tmp_path / "noncanon.fa"
    fp.write_text(">chr1\n" + _wrap("".join(base)))
    return fp


@pytest.fixture
def fasta_reverse_strand(tmp_path: Path) -> Path:
    """600 bp chromosome where the gene lives on the - strand. The
    forward-strand sequence at 101..199 is the reverse complement of the
    canonical ATG..TAA ORF."""
    base = ["A"] * 600
    # Reverse complement of ATG + 32×GCT + TAA = "TTAAGCAGC...AGCAGCCAT"
    from Bio.Seq import Seq
    fwd_orf = "ATG" + "GCT" * 32 + "TAA"   # 102 bp on + strand
    rc = str(Seq(fwd_orf).reverse_complement())
    for i, ch in enumerate(rc):
        base[100 + i] = ch
    fp = tmp_path / "rev.fa"
    fp.write_text(">chr1\n" + _wrap("".join(base)))
    return fp


@pytest.fixture
def gff_at_chrom_edge(tmp_path: Path) -> Path:
    """Gene whose start coordinate is exactly 1 (chromosome edge)."""
    fp = tmp_path / "edge_start.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gE;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=txE;Parent=gE\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=exE;Parent=txE\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tID=cdsE;Parent=txE\n"
    )
    return fp


@pytest.fixture
def gff_one_bp_feature(tmp_path: Path) -> Path:
    """Single-base-pair feature (start == end). NCBI cols 4-5 permit
    start == end as a single-base feature."""
    fp = tmp_path / "single_bp.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t100\t.\t+\t.\tID=g1bp;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t100\t100\t.\t+\t.\tID=tx1bp;Parent=g1bp\n"
        "chr1\ttest\texon\t100\t100\t.\t+\t.\tID=ex1bp;Parent=tx1bp\n"
    )
    return fp


@pytest.fixture
def gff_off_end_feature(tmp_path: Path) -> Path:
    """Feature whose end exceeds typical chromosome length. Used with
    `fasta_short_chrom` (50 bp) to drive the off-end path."""
    fp = tmp_path / "off_end.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t40\t1000\t.\t+\t.\tID=gOff;gene_biotype=protein_coding\n"
        "chr1\ttest\texon\t40\t1000\t.\t+\t.\tID=exOff;Parent=gOff\n"
    )
    return fp


@pytest.fixture
def gff_nested_genes(tmp_path: Path) -> Path:
    """Gene B fully inside an intron of gene A."""
    fp = tmp_path / "nested.gff3"
    fp.write_text(
        "##gff-version 3\n"
        # Outer gene with 2 exons leaving a wide intron 200..300
        "chr1\ttest\tgene\t101\t399\t.\t+\t.\tID=gA;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t399\t.\t+\t.\tID=txA;Parent=gA\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exA1;Parent=txA\n"
        "chr1\ttest\texon\t301\t399\t.\t+\t.\tID=exA2;Parent=txA\n"
        # Inner gene completely inside the intron
        "chr1\ttest\tgene\t220\t280\t.\t-\t.\tID=gB;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t220\t280\t.\t-\t.\tID=txB;Parent=gB\n"
        "chr1\ttest\texon\t220\t280\t.\t-\t.\tID=exB;Parent=txB\n"
    )
    return fp


@pytest.fixture
def gff_shared_exon_isoforms(tmp_path: Path) -> Path:
    """One gene with two mRNA isoforms sharing exon1 via Parent=tx1,tx2."""
    fp = tmp_path / "shared_exon.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t499\t.\t+\t.\tID=gS;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gS\n"
        "chr1\ttest\tmRNA\t101\t499\t.\t+\t.\tID=tx2;Parent=gS\n"
        # Shared first exon
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exShared;Parent=tx1,tx2\n"
        # Distinct downstream exons
        "chr1\ttest\texon\t301\t399\t.\t+\t.\tID=exTx1;Parent=tx1\n"
        "chr1\ttest\texon\t401\t499\t.\t+\t.\tID=exTx2;Parent=tx2\n"
    )
    return fp


@pytest.fixture
def gff_overlapping_opposite_strand(tmp_path: Path) -> Path:
    """Two genes sharing genomic span on opposite strands."""
    fp = tmp_path / "opp_strand.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=gFwd;gene_biotype=protein_coding\n"
        "chr1\ttest\tgene\t150\t350\t.\t-\t.\tID=gRev;gene_biotype=protein_coding\n"
    )
    return fp


@pytest.fixture
def gff_gene_no_transcripts(tmp_path: Path) -> Path:
    """Gene line with no children. Per NCBI annotation data model, this is
    a tolerated state (e.g. gene that is not marked as pseudogene but
    has no child features)."""
    fp = tmp_path / "no_children.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gOnly;gene_biotype=other\n"
    )
    return fp


@pytest.fixture
def gff_prokaryote_gene_to_cds(tmp_path: Path) -> Path:
    """Prokaryote pattern: gene → CDS with no intermediate mRNA, per
    NCBI 'NOTE 2' in the GFF3 specifications."""
    fp = tmp_path / "prok.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=gP;gene_biotype=protein_coding\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cdsP;Parent=gP\n"
    )
    return fp


# ---------------------------------------------------------------------------
# Corruption fixtures (Step 5)
# ---------------------------------------------------------------------------

@pytest.fixture
def gff_missing_directive(tmp_path: Path) -> Path:
    fp = tmp_path / "no_directive.gff3"
    fp.write_text(
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gND;gene_biotype=protein_coding\n"
    )
    return fp


@pytest.fixture
def gff_eight_columns(tmp_path: Path) -> Path:
    """Only 8 tab-separated columns (missing attributes column)."""
    fp = tmp_path / "eight_cols.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\n"
    )
    return fp


@pytest.fixture
def gff_extra_tab_column(tmp_path: Path) -> Path:
    """Ten columns: an extra tab inside the attributes block."""
    fp = tmp_path / "ten_cols.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gEx\textraField=value\n"
    )
    return fp


@pytest.fixture
def gff_crlf_line_endings(tmp_path: Path) -> Path:
    fp = tmp_path / "crlf.gff3"
    fp.write_bytes(
        b"##gff-version 3\r\n"
        b"chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gCRLF;gene_biotype=protein_coding\r\n"
    )
    return fp


@pytest.fixture
def gff_with_bom(tmp_path: Path) -> Path:
    """UTF-8 BOM at file start."""
    fp = tmp_path / "bom.gff3"
    fp.write_bytes(
        b"\xef\xbb\xbf##gff-version 3\n"
        b"chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gBOM;gene_biotype=protein_coding\n"
    )
    return fp


@pytest.fixture
def gff_dangling_parent(tmp_path: Path) -> Path:
    """mRNA references gene 'ghost' that does not appear in the file."""
    fp = tmp_path / "dangling.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tmRNA\t100\t200\t.\t+\t.\tID=txD;Parent=ghost\n"
    )
    return fp


@pytest.fixture
def gff_unencoded_semicolon(tmp_path: Path) -> Path:
    """Reserved character ';' inside an attribute value, NOT
    percent-encoded. Per NCBI, this should be %3B."""
    fp = tmp_path / "unencoded.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gU;Note=alpha;beta\n"
    )
    return fp


@pytest.fixture
def gff_negative_coords(tmp_path: Path) -> Path:
    fp = tmp_path / "neg.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t-5\t100\t.\t+\t.\tID=gNeg\n"
    )
    return fp


@pytest.fixture
def gff_truncated_no_newline(tmp_path: Path) -> Path:
    fp = tmp_path / "trunc.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gT;gene_biotype=protein_coding"
    )  # no trailing newline
    return fp


@pytest.fixture
def gff_duplicate_id_collision(tmp_path: Path) -> Path:
    """Same ID on two rows whose (seqid, type) differ — true collision per
    NCBI ID semantics (multi-row same-feature is only legal when all
    fields match)."""
    fp = tmp_path / "dup_id.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=dup;gene_biotype=protein_coding\n"
        "chr2\ttest\tmRNA\t300\t400\t.\t+\t.\tID=dup;Parent=somegene\n"
    )
    return fp


@pytest.fixture(autouse=True)
def _scrub_gffutils_db_cache(tmp_path):
    """gffutils writes <gff>_db next to the GFF; tmp_path is per-test so
    no manual cleanup is normally needed, but defensively wipe any
    leftovers from prior failed runs in the repo root."""
    yield
    for stale in REPO_ROOT.glob("**/*.gff3_db"):
        # Only delete inside tmp / tests, never in source data
        if "tests" in stale.parts and stale.is_file():
            try:
                stale.unlink()
            except OSError:
                pass
