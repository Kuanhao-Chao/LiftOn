"""Unit tests for Phase D — duplicate-safe cross-locus miniprot rescue
(``lifton.cross_locus_rescue``). These prove the env-independent CORRECTNESS:
the score index, the duplicate-safe GFF3 rewrite, the weak/cross-locus/min-gain
gating, and flag-off inertness. The aggregate biological win is a separate A/B.
"""
from types import SimpleNamespace

import pytest

from lifton import cross_locus_rescue as xlocus, lifton_class


# ── a minimal valid gene block on chrB used as the miniprot replacement text ──
REPLACEMENT_TEXT = (
    "chrB\tLiftOn\tgene\t100\t400\t.\t+\t.\tID=gene-X;Name=X\n"
    "chrB\tLiftOn\tmRNA\t100\t400\t.\t+\t.\tID=rna-X;Parent=gene-X;lifton_rescue=cross_locus\n"
    "chrB\tLiftOn\texon\t100\t400\t.\t+\t.\tID=exon-X-1;Parent=rna-X\n"
    "chrB\tLiftOn\tCDS\t100\t400\t.\t+\t0\tID=cds-X;Parent=rna-X\n"
)

OUTPUT_GFF3 = (
    "##gff-version 3\n"
    "#!some-directive\n"
    # gene-X (weak, on chrA) — should be DROPPED on replacement
    "chrA\tLiftOn\tgene\t50\t90\t.\t-\t.\tID=gene-X;Name=X\n"
    "chrA\tLiftOn\tmRNA\t50\t90\t.\t-\t.\tID=rna-X;Parent=gene-X\n"
    "chrA\tLiftOn\texon\t50\t90\t.\t-\t.\tID=exon-X-1;Parent=rna-X\n"
    "chrA\tLiftOn\tCDS\t50\t90\t.\t-\t0\tID=cds-X;Parent=rna-X\n"
    # gene-Y (strong) — must be UNTOUCHED
    "chrA\tLiftOn\tgene\t500\t900\t.\t+\t.\tID=gene-Y;Name=Y\n"
    "chrA\tLiftOn\tmRNA\t500\t900\t.\t+\t.\tID=rna-Y;Parent=gene-Y\n"
    "chrA\tLiftOn\texon\t500\t900\t.\t+\t.\tID=exon-Y-1;Parent=rna-Y\n"
    "chrA\tLiftOn\tCDS\t500\t900\t.\t+\t0\tID=cds-Y;Parent=rna-Y\n"
)


# ───────────────────────────── small units ───────────────────────────────────
def test_index_emitted_keys_by_ref_gene(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text(
        # tid liftoff miniprot dna AA annotation status loc
        "rna-X\t0\t0\t0\t0.10\tLiftoff\t-\tchrA:50-90\n"
        "rna-X_1\t0\t0\t0\t0.20\tLiftoff\t-\tchrC:5-9\n"     # a copy -> base rna-X
        "rna-Y\t0\t0\t0\t0.95\tLiftoff\t-\tchrA:500-900\n"
    )
    rev = {"rna-X": "gene-X", "rna-Y": "gene-Y"}
    idx = xlocus._index_emitted(str(score), {}, rev)
    assert idx["gene-X"]["best_pi"] == pytest.approx(0.20)        # max over copies
    assert set(idx["gene-X"]["intervals"]) == {("chrA", 50, 90), ("chrC", 5, 9)}
    assert idx["gene-Y"]["best_pi"] == pytest.approx(0.95)


def test_overlaps_any_requires_same_seqid_and_coordinate_overlap():
    # Two different genes on the SAME chromosome at DIFFERENT positions are NOT
    # "the same locus" -- a same-seqid-only check would wrongly treat any
    # same-chromosome paralog as already-covered (the bug the integration test
    # caught: a synthetic same-chromosome cross-locus candidate was rejected).
    ivs = [("chrA", 50, 90)]
    assert xlocus._overlaps_any("chrA", 50, 90, ivs) is True          # exact overlap
    assert xlocus._overlaps_any("chrA", 60, 70, ivs) is True          # nested overlap
    assert xlocus._overlaps_any("chrA", 200, 300, ivs) is False       # same seqid, no overlap
    assert xlocus._overlaps_any("chrB", 50, 90, ivs) is False         # different seqid


def test_rewrite_drops_replaced_block_and_appends(tmp_path):
    out = tmp_path / "out.gff3"
    out.write_text(OUTPUT_GFF3)
    n_blocks, n_lines = xlocus._rewrite_output(
        str(out), {"gene-X"}, {"gene-X": REPLACEMENT_TEXT}, {})
    text = out.read_text()
    assert n_blocks == 1 and n_lines == 4              # one 4-line gene-X block dropped
    # gene-X appears EXACTLY once now (the appended replacement) -> duplicate-safe
    assert text.count("ID=gene-X;") == 1
    assert "chrB\tLiftOn\tgene" in text               # replacement is on chrB
    assert "lifton_rescue=cross_locus" in text
    # gene-Y untouched + directives preserved
    assert text.count("ID=gene-Y;") == 1
    assert text.startswith("##gff-version 3\n#!some-directive\n")


def test_resolution_env_and_args(monkeypatch):
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    assert xlocus._enabled(SimpleNamespace(cross_locus_rescue=False)) is False
    assert xlocus._enabled(SimpleNamespace(cross_locus_rescue=True)) is True
    monkeypatch.setenv("LIFTON_CROSS_LOCUS_RESCUE", "1")
    assert xlocus._enabled(SimpleNamespace(cross_locus_rescue=False)) is True
    monkeypatch.setenv("LIFTON_CROSS_LOCUS_RESCUE", "0")
    assert xlocus._enabled(SimpleNamespace(cross_locus_rescue=True)) is False
    monkeypatch.setenv("LIFTON_CROSS_LOCUS_MAX_LIFTOFF", "0.4")
    monkeypatch.setenv("LIFTON_CROSS_LOCUS_MIN_GAIN", "0.25")
    assert xlocus._max_liftoff(SimpleNamespace()) == pytest.approx(0.4)
    assert xlocus._min_gain(SimpleNamespace()) == pytest.approx(0.25)


# ─────────────────── full pass (mocked model builder) ─────────────────────────
class _FakeMtrans:
    def __init__(self, tid, seqid, start, end):
        self.attributes = {"ID": [tid]}
        self.seqid = seqid
        self.start = start
        self.end = end
        self.id = tid


class _FakeMDB:
    def __init__(self, mtranscripts):
        self._m = mtranscripts

    def features_of_type(self, ftype):
        return list(self._m)

    def children(self, mtrans, featuretype=None):
        return []          # multi-exon CDS check -> not a processed pseudogene


class _FakeReplGene:
    def __init__(self, tid):
        self.entry = SimpleNamespace(id="gene-X", attributes={})
        self.transcripts = {tid: SimpleNamespace(
            entry=SimpleNamespace(id=tid, attributes={}))}

    def orf_search_protein(self, tid, ref_tid, tgt_fai, ref_proteins, ref_trans, status):
        status.lifton_aa = 0.9
        return ("", 0.9)

    def add_lifton_gene_status_attrs(self, source):
        pass

    def add_lifton_trans_status_attrs(self, tid, status):
        pass

    def write_entry(self, fw, stats):
        fw.write(REPLACEMENT_TEXT)


def _patch_builder(monkeypatch, pi_gene_id="gene-X", mp_seqid="chrB", ref_gene="gene-X"):
    monkeypatch.setattr(
        xlocus.lifton_utils, "get_ref_ids_miniprot",
        lambda rev, mid, m2r: (ref_gene, "rna-X"))

    def _fake_build(mtrans, mdb, ref_db, rg, rt, tgt_fai, refp, reft, tree, rfd, args):
        return (_FakeReplGene("rna-X"),
                SimpleNamespace(entry=SimpleNamespace(id="rna-X")),
                "rna-X", lifton_class.Lifton_Status())
    monkeypatch.setattr(xlocus.run_miniprot, "lifton_miniprot_with_ref_protein",
                        _fake_build)


def _run_pass(tmp_path, monkeypatch, *, mp_seqid="chrB", mp_start=100, mp_end=400,
              flag=True, max_liftoff=0.5, min_gain=0.3):
    out = tmp_path / "out.gff3"; out.write_text(OUTPUT_GFF3)
    score = tmp_path / "score.txt"
    score.write_text(
        "rna-X\t0\t0\t0\t0.10\tLiftoff\t-\tchrA:50-90\n"
        "rna-Y\t0\t0\t0\t0.95\tLiftoff\t-\tchrA:500-900\n")
    _patch_builder(monkeypatch)
    mtrans = _FakeMtrans("MP1", mp_seqid, mp_start, mp_end)
    args = SimpleNamespace(cross_locus_rescue=flag,
                           cross_locus_max_liftoff=max_liftoff,
                           cross_locus_min_gain=min_gain, debug=False)
    rev = {"rna-X": "gene-X", "rna-Y": "gene-Y"}
    n = xlocus.cross_locus_rescue_pass(
        str(out), str(score), _FakeMDB([mtrans]),
        SimpleNamespace(db_connection=None), {}, None, {"rna-X": "P"},
        {"rna-X": "T"}, {}, {"MP1": "rna-X"}, {}, {}, rev, args)
    return n, out.read_text()


def test_full_pass_replaces_weak_cross_locus_gene(tmp_path, monkeypatch):
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MAX_LIFTOFF", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MIN_GAIN", raising=False)
    n, text = _run_pass(tmp_path, monkeypatch)
    assert n == 1
    assert text.count("ID=gene-X;") == 1               # no duplicate
    assert "chrB\tLiftOn\tgene" in text                # moved to miniprot locus
    assert text.count("ID=gene-Y;") == 1               # strong gene untouched


def test_full_pass_skips_same_locus(tmp_path, monkeypatch):
    # miniprot genuinely OVERLAPS the weak lift's interval (chrA:50-90) ->
    # candidate-3's job, not this.
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MAX_LIFTOFF", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MIN_GAIN", raising=False)
    n, text = _run_pass(tmp_path, monkeypatch, mp_seqid="chrA", mp_start=60, mp_end=120)
    assert n == 0
    assert text == OUTPUT_GFF3                          # untouched


def test_full_pass_rescues_same_seqid_different_position(tmp_path, monkeypatch):
    # REGRESSION: two genes on the SAME chromosome at DIFFERENT (non-overlapping)
    # positions are NOT "the same locus" -- an earlier version of this gate
    # checked raw seqid membership instead of coordinate overlap, so it wrongly
    # rejected this exact case (caught by
    # tests/test_integration_pipeline.py::TestCrossLocusRescuePass, whose
    # synthetic single-chromosome fixture exercises precisely this).
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MAX_LIFTOFF", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MIN_GAIN", raising=False)
    # chrA:1000-1300 -- same chromosome as gene-X's weak lift (chrA:50-90), but
    # clear of both that interval AND unrelated gene-Y's (chrA:500-900).
    n, text = _run_pass(tmp_path, monkeypatch, mp_seqid="chrA", mp_start=1000, mp_end=1300)
    assert n == 1
    assert text.count("ID=gene-X;") == 1
    assert "chrB\tLiftOn\tgene" in text


def test_full_pass_skips_strong_gene(tmp_path, monkeypatch):
    # raise the weak threshold target so even the 0.10 gene counts as 'strong'
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MAX_LIFTOFF", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MIN_GAIN", raising=False)
    n, text = _run_pass(tmp_path, monkeypatch, max_liftoff=0.05)
    assert n == 0
    assert text == OUTPUT_GFF3


def test_full_pass_respects_min_gain(tmp_path, monkeypatch):
    # require a gain so large the 0.10 -> 0.90 jump (0.80) does not clear it
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MAX_LIFTOFF", raising=False)
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_MIN_GAIN", raising=False)
    n, text = _run_pass(tmp_path, monkeypatch, min_gain=0.95)
    assert n == 0
    assert text == OUTPUT_GFF3


def test_flag_off_is_inert(tmp_path, monkeypatch):
    monkeypatch.delenv("LIFTON_CROSS_LOCUS_RESCUE", raising=False)
    n, text = _run_pass(tmp_path, monkeypatch, flag=False)
    assert n == 0
    assert text == OUTPUT_GFF3                          # byte-identical when OFF
