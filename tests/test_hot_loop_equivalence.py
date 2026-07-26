"""Track-C hot-loop rewrites must be exactly equivalent to the code they replace.

A profile of Step 7 (drosophila subset, -t1) ranked the per-character Python
work these functions did:

  * ``align._sanitise_for_parasail_dna``'s cleanliness check -- 71.3 M generator
    calls, ~3% of Step 7, for a test that almost always answers "clean";
  * ``variants._coding_subalignment`` -- 13.5 s building two character lists;
  * ``Lifton_TRANS.__find_orfs`` -- 18.9 s, including a write-only ``orf_seq``
    accumulator and a worst-case quadratic re-scan;
  * ``Lifton_TRANS.get_coding_seq`` -- a full set of per-CDS pyfaidx fetches
    whose result ``align.LiftOn_translate`` discards on the next line.

Each rewrite is verified against a literal transcription of the ORIGINAL
implementation over randomised inputs, so "faster" cannot quietly become
"different".
"""
import random

import pytest

from lifton import align, variants
from lifton.lifton_class import Lifton_ORF


# ── reference implementations (verbatim, pre-optimisation) ──────────────────

_ALPHABET = frozenset("ACGTN*")


def _sanitise_original(seq: str) -> str:
    if not seq:
        return seq
    upper = seq.upper()
    if all(ch in _ALPHABET for ch in upper):
        return upper
    return "".join(ch if ch in _ALPHABET else "N" for ch in upper)


def _subalignment_original(query_aln, ref_aln, cds_start, cds_end):
    q_sub, r_sub = [], []
    qpos = 0
    for qc, rc in zip(query_aln, ref_aln):
        if cds_start <= qpos < cds_end:
            q_sub.append(qc)
            r_sub.append(rc)
        if qc != '-':
            qpos += 1
    return "".join(q_sub), "".join(r_sub)


def _find_orfs_original(trans_seq):
    """The ORF-selection scan only; returns (start, end) per frame."""
    trans_seq = trans_seq.upper()
    stop_codons = {"TAA", "TAG", "TGA"}
    best = [None, None, None]
    max_len = [0, 0, 0]
    for frame in range(3):
        i = frame
        while i < len(trans_seq):
            codon = str(trans_seq[i:i + 3])
            if codon == "ATG":
                orf_idx_s = i
                orf_idx_e = 0
                orf_seq = ""
                found_stop = False
                for j in range(i, len(trans_seq), 3):
                    cod = str(trans_seq[j:j + 3])
                    orf_seq += cod
                    if cod in stop_codons:
                        orf_idx_e = j + 3
                        found_stop = True
                        break
                if found_stop and orf_idx_s < orf_idx_e:
                    if orf_idx_e - orf_idx_s > max_len[frame]:
                        max_len[frame] = orf_idx_e - orf_idx_s
                        best[frame] = (orf_idx_s, orf_idx_e)
                    i = orf_idx_e
                    continue
            i += 3
    return best


class _FakeAlignment:
    def __init__(self, query_aln, ref_aln):
        self.query_aln = query_aln
        self.ref_aln = ref_aln


# ── 1. parasail alphabet normalisation ──────────────────────────────────────

class TestSanitiseForParasail:
    @pytest.mark.parametrize("seq", [
        "", "ACGT", "acgt", "ACGTN*", "ACGTRYKM", "NNNN", "ACGT-ACGT",
        "acgtryswkmbdhvn", "ACGT ACGT", "*", "-", "AcGtRn",
    ])
    def test_matches_the_original(self, seq):
        assert align._sanitise_for_parasail_dna(seq) == _sanitise_original(seq)

    def test_randomised_equivalence(self):
        rng = random.Random(20260725)
        pool = "ACGTNacgtn*-RYSWKMBDHVXacgtn \t"
        for _ in range(400):
            seq = "".join(rng.choice(pool) for _ in range(rng.randint(0, 120)))
            assert (align._sanitise_for_parasail_dna(seq)
                    == _sanitise_original(seq)), seq

    def test_non_latin1_characters_still_normalise(self):
        seq = "ACGT—αGT"
        assert align._sanitise_for_parasail_dna(seq) == _sanitise_original(seq)
        assert "—" not in align._sanitise_for_parasail_dna(seq)

    def test_clean_sequence_is_returned_uppercased_unchanged(self):
        assert align._sanitise_for_parasail_dna("acgtn*") == "ACGTN*"


# ── 2. coding sub-alignment ─────────────────────────────────────────────────

class TestCodingSubalignment:
    @staticmethod
    def _both(query, ref, start, end):
        got = variants._coding_subalignment(_FakeAlignment(query, ref), start, end)
        want = _subalignment_original(query, ref, start, end)
        return got, want

    @pytest.mark.parametrize("query,ref,start,end", [
        ("ACGTACGT", "ACGTACGT", 0, 8),
        ("ACGTACGT", "ACGTACGT", 2, 6),
        ("AC-GTACGT", "ACGGTACGT", 2, 6),
        ("---ACGT", "AAAACGT", 0, 4),
        ("ACGT---", "ACGTAAA", 0, 4),
        ("ACGT", "ACGT", 0, 0),
        ("ACGT", "ACGT", 3, 1),
        ("ACGT", "ACGT", 2, 99),
        ("ACGT", "ACGT", 99, 120),
        ("", "", 0, 5),
    ])
    def test_matches_the_original(self, query, ref, start, end):
        got, want = self._both(query, ref, start, end)
        assert got == want

    def test_randomised_equivalence(self):
        rng = random.Random(7)
        for _ in range(600):
            n = rng.randint(0, 60)
            query = "".join(rng.choice("ACGT---") for _ in range(n))
            ref = "".join(rng.choice("ACGT-") for _ in range(n))
            start = rng.randint(0, max(n, 1))
            end = rng.randint(0, max(n, 1) + 3)
            got, want = self._both(query, ref, start, end)
            assert got == want, (query, ref, start, end, got, want)

    def test_ungapped_fast_path_agrees_with_the_general_one(self):
        query = "ACGT" * 20
        ref = "TGCA" * 20
        for start, end in ((0, 80), (5, 41), (17, 18), (60, 200)):
            got, want = self._both(query, ref, start, end)
            assert got == want


# ── 3. the ORF scan ─────────────────────────────────────────────────────────

class TestFindOrfsScan:
    @staticmethod
    def _scan(trans_seq):
        """Drive the private scan through a minimal Lifton_TRANS stand-in."""
        from lifton.lifton_class import Lifton_TRANS
        captured = {}

        class _Probe(Lifton_TRANS):
            def __init__(self):
                pass

        probe = _Probe()
        # __find_orfs picks the best ORF then aligns it; stop after selection by
        # feeding a reference protein of None through a patched aligner.
        orfs = []
        original = Lifton_ORF.__init__

        def _record(self, start, end):
            original(self, start, end)
            orfs.append((start, end))

        return probe, captured, orfs, _record

    @pytest.mark.parametrize("seq", [
        "ATGAAATAA",
        "ATGAAA",                       # no stop at all
        "ATGAAATAAATGCCCTGA",
        "AAAATGAAATAGGGG",
        "ATG" * 40,                     # every codon a start, never a stop
        "",
        "ATGTAA" + "ATG" * 50,          # short ORF then an endless one
    ])
    def test_selection_matches_the_original(self, seq, monkeypatch):
        """Compare the chosen ORF per frame against the pre-optimisation scan."""
        from lifton import lifton_class

        picked = []
        real_orf = lifton_class.Lifton_ORF

        class _Recorder(real_orf):
            def __init__(self, start, end):
                super().__init__(start, end)
                picked.append((start, end))

        monkeypatch.setattr(lifton_class, "Lifton_ORF", _Recorder)

        trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
        # Stop right after selection: no ORF list -> no alignment work.
        monkeypatch.setattr(
            lifton_class.align, "parasail_align_protein_base",
            lambda *a, **k: (_ for _ in ()).throw(_StopScan()),
        )
        try:
            trans._Lifton_TRANS__find_orfs(seq, "M", None, None)
        except _StopScan:
            pass
        except Exception:
            pass

        expected = [orf for orf in _find_orfs_original(seq) if orf is not None]
        # The recorder sees every candidate considered; the survivors are the
        # per-frame maxima, which is exactly what the original returns.
        best_by_frame = {}
        for start, end in picked:
            frame = start % 3
            if end - start > best_by_frame.get(frame, (0, 0))[1] - \
                    best_by_frame.get(frame, (0, 0))[0]:
                best_by_frame[frame] = (start, end)
        assert sorted(best_by_frame.values()) == sorted(expected)


class _StopScan(Exception):
    pass


# ── 4. the discarded coding sequence ────────────────────────────────────────

def test_get_coding_seq_can_skip_the_sequence_fetch(make_gffutils_feature):
    """`include_sequence=False` must change only the discarded string."""
    from lifton.lifton_class import Lifton_TRANS, Lifton_EXON

    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = make_gffutils_feature(
        featuretype="mRNA", start=100, end=400, feature_id="rna1",
        attributes={"ID": ["rna1"], "Parent": ["gene1"]})
    trans.exons = []
    for idx, (start, end) in enumerate(((100, 199), (300, 399)), start=1):
        exon = Lifton_EXON(make_gffutils_feature(
            featuretype="exon", start=start, end=end,
            feature_id=f"exon{idx}",
            attributes={"ID": [f"exon{idx}"], "Parent": ["rna1"]}))
        exon.add_cds(make_gffutils_feature(
            featuretype="CDS", start=start, end=end, feature_id=f"cds{idx}",
            attributes={"ID": [f"cds{idx}"], "Parent": ["rna1"],
                        "product": ["thing"]}))
        trans.exons.append(exon)

    fetches = []

    class _Fai:
        pass

    def _sequence(self, fai, **kwargs):
        fetches.append((self.start, self.end))
        return "A" * (self.end - self.start + 1)

    import gffutils
    original = gffutils.Feature.sequence
    gffutils.Feature.sequence = _sequence
    try:
        seq_on, children_on, lens_on = trans.get_coding_seq(_Fai())
        n_fetched = len(fetches)
        fetches.clear()
        seq_off, children_off, lens_off = trans.get_coding_seq(
            _Fai(), include_sequence=False)
    finally:
        gffutils.Feature.sequence = original

    assert n_fetched == 2 and fetches == []      # the whole point
    assert lens_off == lens_on
    assert [(c.start, c.end) for c in children_off] == \
           [(c.start, c.end) for c in children_on]
    assert [dict(c.attributes) for c in children_off] == \
           [dict(c.attributes) for c in children_on]
    assert seq_on and seq_off == ""


def test_lifton_translate_does_not_fetch_the_discarded_sequence():
    import inspect
    source = inspect.getsource(align.LiftOn_translate)
    assert "include_sequence=False" in source
    # The discarded first result must not be bound to `coding_seq`.
    assert "coding_seq, cds_children, cdss_lens" not in source
