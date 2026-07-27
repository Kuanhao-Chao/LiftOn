"""Unit tests for candidate-3 — the standalone miniprot-only best-of-outcome
merge candidate in ``run_liftoff.process_liftoff_with_protein`` (Figure-4 fix).

The 2-way best-of-outcome merge keeps max(Liftoff+ORF, chained-merge+ORF) but
never considers miniprot's NATIVE model; at very-distant divergence that 2-way
winner can collapse to a truncated Liftoff stub while miniprot's clean model is
discarded. Candidate-3 evaluates miniprot's native CDS + ORF-rescue and emits it
ONLY when STRICTLY better, so it is additive (per-transcript identity
non-decreasing) and a no-op wherever miniprot is not strictly better.

These drive the 3-way decision directly with controlled candidate outcomes:
``Lifton_GENE.orf_search_protein`` is replaced by a stateful queue that yields a
chosen identity per call (call 1 = merge, call 2 = Liftoff, call 3 = miniprot),
and the two alignment helpers + the chaining are stubbed. This exercises the
real snapshot/restore + tie-break logic in ``process_liftoff_with_protein``
without needing a synthetic genome where the merge collapses just so.
"""
from types import SimpleNamespace

import pytest

from lifton import run_liftoff, lifton_class


# ─── lightweight fakes (only the attributes the decision logic reads) ──────────
class _FakeEntry:
    def __init__(self, fid, start, end, ftype="mRNA", strand="+"):
        self.id = fid
        self.start = start
        self.end = end
        self.strand = strand
        self.source = "Liftoff"
        self.featuretype = ftype
        self.attributes = {}


class _FakeAln:
    """Stands in for a Lifton alignment — only ``identity``/``cds_children`` used."""
    def __init__(self, identity, cds_children):
        self.identity = identity
        self.cds_children = cds_children


class _FakeRawCDS:
    """A raw gffutils-Feature-shaped CDS that ``Lifton_CDS(c)`` can wrap."""
    def __init__(self, start, end, strand="+"):
        self.start = start
        self.end = end
        self.strand = strand
        self.source = "miniprot"
        self.featuretype = "CDS"
        self.attributes = {}


class _FakeCDSWrap:
    """A chained-CDS entry — needs ``.entry.start/.end`` for the _merge_noop test."""
    def __init__(self, start, end):
        self.entry = _FakeEntry("cds", start, end, ftype="CDS")


class _FakeGene:
    def __init__(self, orf_queue):
        self.entry = _FakeEntry("gene1", 1, 1000, ftype="gene")
        self.orf_queue = list(orf_queue)
        self.orf_calls = 0
        self.update_calls = []          # list of [(start,end), ...] applied
        self.added_exons = []           # clean-rebuild: (start,end) per add_exon call
        self.added_cds = []             # clean-rebuild: (start,end) per add_cds call
        self.boundary_calls = 0

    def update_cds_list(self, tid, cds_list, optimize=False):
        self.update_calls.append([(c.entry.start, c.entry.end) for c in cds_list])

    def add_exon(self, tid, c):
        self.added_exons.append((c.start, c.end))

    def add_cds(self, tid, c):
        self.added_cds.append((c.start, c.end))

    def update_boundaries(self):
        self.boundary_calls += 1

    def orf_search_protein(self, tid, ref_tid, tgt_fai, ref_proteins, ref_trans, status):
        v = self.orf_queue[self.orf_calls]
        self.orf_calls += 1
        status.lifton_aa = v
        return ("", v)


class _FakeTrans:
    def __init__(self):
        self.entry = _FakeEntry("mrna1", 10, 900)
        self.exons = []

    def update_boundaries(self):
        pass


def _drive(monkeypatch, orf_queue, *, liftoff_id=0.5,
           liftoff_cds=(50, 140), chain_cds=(50, 150), mini_cds=(200, 300),
           mini_strand="+"):
    """Run ``process_liftoff_with_protein`` with stubbed alignments + a stateful
    ORF queue. Returns (lifton_status, gene). ``mini_strand`` sets the strand of
    the miniprot CDS children (the transcript is always ``+``), so passing ``-``
    exercises the opposite-strand guard."""
    gene = _FakeGene(orf_queue)
    trans = _FakeTrans()
    status = lifton_class.Lifton_Status()

    monkeypatch.setattr(
        run_liftoff.lifton_utils, "LiftOn_liftoff_alignment",
        lambda *a, **k: _FakeAln(liftoff_id, [_FakeRawCDS(*liftoff_cds)]))
    monkeypatch.setattr(
        run_liftoff.lifton_utils, "LiftOn_miniprot_alignment",
        lambda *a, **k: (_FakeAln(0.9, [_FakeRawCDS(*mini_cds, strand=mini_strand)]), True))
    monkeypatch.setattr(
        run_liftoff.protein_maximization, "chaining_algorithm",
        lambda *a, **k: ([_FakeCDSWrap(*chain_cds)], []))

    locus = SimpleNamespace(seqid="chr1")
    run_liftoff.process_liftoff_with_protein(
        locus, gene, trans, {}, None, {}, None, "ref1", {}, {},
        None, False, status, False)
    return status, gene


def test_candidate3_emits_miniprot_when_strictly_better(monkeypatch):
    # merge=0.5, Liftoff=0.5 (2-way winner=merge), miniprot=0.9 → miniprot wins.
    status, gene = _drive(monkeypatch, [0.5, 0.5, 0.9])
    assert status.annotation == "LiftOn_miniprot"
    assert status.lifton_aa == pytest.approx(0.9)
    assert gene.orf_calls == 3                       # merge, Liftoff, miniprot
    # candidate-3 rebuilds a CLEAN scaffold (add_exon/add_cds), NOT a reconciling
    # update_cds_list graft (the corruption the chicken A/B caught).
    assert gene.added_exons == [(200, 300)]
    assert gene.added_cds == [(200, 300)]
    assert gene.boundary_calls == 1


def test_candidate3_reverts_to_2way_when_not_better(monkeypatch):
    # miniprot ties the 2-way winner (0.5) → NOT strictly better → revert.
    status, gene = _drive(monkeypatch, [0.5, 0.5, 0.5])
    assert status.annotation == "LiftOn_chaining_algorithm"
    assert status.lifton_aa == pytest.approx(0.5)
    assert gene.orf_calls == 3                       # candidate-3 ran but reverted


def test_candidate3_disabled_by_env_restores_2way(monkeypatch):
    monkeypatch.setenv("LIFTON_MINIPROT_CANDIDATE", "0")
    status, gene = _drive(monkeypatch, [0.5, 0.5])   # only 2 ORF calls expected
    assert status.annotation == "LiftOn_chaining_algorithm"
    assert status.lifton_aa == pytest.approx(0.5)
    assert gene.orf_calls == 2                       # candidate-3 never ran


def test_candidate3_skipped_when_merge_already_perfect(monkeypatch):
    # merge reaches 1.0 → optimize-fast keeps it; _best_outcome>=1.0 → candidate-3
    # guard (best<1.0) is False → skipped, single ORF call.
    status, gene = _drive(monkeypatch, [1.0])
    assert status.annotation == "LiftOn_chaining_algorithm"
    assert status.lifton_aa == pytest.approx(1.0)
    assert gene.orf_calls == 1


def test_candidate3_skipped_when_no_miniprot_children(monkeypatch):
    # miniprot has no cds_children → candidate-3 cannot fire (the cross-locus
    # case): 2-way winner stands.
    gene = _FakeGene([0.5, 0.5])
    trans = _FakeTrans()
    status = lifton_class.Lifton_Status()
    monkeypatch.setattr(
        run_liftoff.lifton_utils, "LiftOn_liftoff_alignment",
        lambda *a, **k: _FakeAln(0.5, [_FakeRawCDS(50, 140)]))
    monkeypatch.setattr(
        run_liftoff.lifton_utils, "LiftOn_miniprot_alignment",
        lambda *a, **k: (_FakeAln(0.9, []), True))      # empty cds_children
    monkeypatch.setattr(
        run_liftoff.protein_maximization, "chaining_algorithm",
        lambda *a, **k: ([_FakeCDSWrap(50, 150)], []))
    run_liftoff.process_liftoff_with_protein(
        SimpleNamespace(seqid="chr1"), gene, trans, {}, None, {}, None,
        "ref1", {}, {}, None, False, status, False)
    assert status.annotation == "LiftOn_chaining_algorithm"
    assert gene.orf_calls == 2


def test_merge_and_candidate3_skipped_when_miniprot_opposite_strand(monkeypatch):
    # GitHub issue #33: a miniprot hit on the OPPOSITE strand from the lifted
    # transcript (transcript "+", miniprot CDS "-") is an antisense/spurious match.
    # Chaining its CDS into the merge under the "+" mRNA emits strand-inconsistent
    # GFF3 (exon "+" / CDS "-"; the drosophila rna-NM_176527.1 failure: 0->4
    # strand_consistency errors, neutral protein id 0.138->0.045). The strand guard
    # now skips the ENTIRE miniprot merge — not just candidate-3 (which the pre-#33
    # guard did): process_liftoff_with_protein makes NO ORF call and leaves orf_done
    # False, so the caller falls back to pure Liftoff + ORF-rescue for this transcript.
    status, gene = _drive(monkeypatch, [0.5, 0.5], mini_strand="-")
    assert status.annotation != "LiftOn_chaining_algorithm"   # antisense merge did NOT run
    assert gene.orf_calls == 0                                 # merge block entirely skipped
    assert gene.added_exons == []                             # no miniprot model built


def test_candidate3_fires_when_miniprot_same_strand(monkeypatch):
    # Control for the guard: same strand ("+"/"+") + strictly-better miniprot
    # (0.9 > 0.5) still fires candidate-3, exactly as without the guard.
    status, gene = _drive(monkeypatch, [0.5, 0.5, 0.9], mini_strand="+")
    assert status.annotation == "LiftOn_miniprot"
    assert gene.orf_calls == 3
    assert gene.added_exons == [(200, 300)]


# ═══════════════════════════════════════════════════════════════════════════
# Regression test with REAL Lifton_GENE/Lifton_TRANS/Lifton_EXON/Lifton_CDS
# machinery (not the mocked fakes above). Proves the clean-rebuild fix: when
# miniprot's CDS boundaries only PARTIALLY overlap the pre-existing Liftoff
# exon scaffold (the exact non-nesting condition that made the old
# update_cds_list-based graft corrupt the structure — stale + duplicate exons,
# caught by the chicken A/B: 72 catastrophic regressions up to 0.809 -> 0.004),
# the rebuilt transcript contains EXACTLY miniprot's CDS-only exons: no stale
# Liftoff-only exon, no duplicates. Only Lifton_GENE.orf_search_protein is
# stubbed (a controllable identity queue); every exon/CDS mutation — including
# candidate 2's real update_cds_list reconciliation — runs for real.
# ═══════════════════════════════════════════════════════════════════════════
class TestCandidate3CleanRebuildRealObjects:
    def test_no_stale_or_duplicate_exons_after_clean_rebuild(
            self, make_gffutils_feature, fake_args, ref_features_dict_one_gene,
            monkeypatch):
        gene_entry = make_gffutils_feature(
            featuretype="gene", start=101, end=900, feature_id="gene1")
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1", gffutil_entry_gene=gene_entry,
            ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["protein_coding"]},
            tree_dict={}, ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
        )
        trans_entry = make_gffutils_feature(
            featuretype="mRNA", start=101, end=900, feature_id="tx1")
        trans = gene.add_transcript(
            "tx1", trans_entry, {"ID": ["tx1"], "Parent": ["gene1"]})
        # Pre-existing (poor) Liftoff exon/CDS scaffold.
        trans.add_exon(make_gffutils_feature(
            featuretype="exon", start=101, end=199, feature_id="exon1"))
        trans.add_exon(make_gffutils_feature(
            featuretype="exon", start=301, end=399, feature_id="exon2"))
        trans.add_cds(make_gffutils_feature(
            featuretype="CDS", start=101, end=199, feature_id="cds1"))
        trans.add_cds(make_gffutils_feature(
            featuretype="CDS", start=301, end=399, feature_id="cds2"))
        assert len(trans.exons) == 2  # sanity: the pre-existing scaffold

        # miniprot's CDS PARTIALLY overlaps the Liftoff exons (150-250 overlaps
        # exon1's tail; 350-450 overlaps exon2's tail) -- the non-nesting
        # condition that made the old graft leave stale/duplicate exons behind.
        mini_cds_1 = make_gffutils_feature(
            featuretype="CDS", start=150, end=250, feature_id="mp_cds1")
        mini_cds_2 = make_gffutils_feature(
            featuretype="CDS", start=350, end=450, feature_id="mp_cds2")

        monkeypatch.setattr(
            run_liftoff.lifton_utils, "LiftOn_liftoff_alignment",
            lambda *a, **k: _FakeAln(0.5, [
                make_gffutils_feature(featuretype="CDS", start=101, end=199),
                make_gffutils_feature(featuretype="CDS", start=301, end=399),
            ]))
        monkeypatch.setattr(
            run_liftoff.lifton_utils, "LiftOn_miniprot_alignment",
            lambda *a, **k: (_FakeAln(0.9, [mini_cds_1, mini_cds_2]), True))
        monkeypatch.setattr(
            run_liftoff.protein_maximization, "chaining_algorithm",
            # Deliberately shift the 2nd segment's end by 1 vs liftoff_aln's
            # cds_children so _merge_noop is False and all 3 candidates run.
            lambda *a, **k: ([
                lifton_class.Lifton_CDS(make_gffutils_feature(
                    featuretype="CDS", start=101, end=199)),
                lifton_class.Lifton_CDS(make_gffutils_feature(
                    featuretype="CDS", start=301, end=398)),
            ], []))

        orf_queue = [0.5, 0.5, 0.9]   # merge, liftoff, miniprot -> miniprot wins
        calls = {"n": 0}

        def _stub_orf(self, tid, ref_tid, tgt_fai, ref_proteins, ref_trans, status):
            v = orf_queue[calls["n"]]
            calls["n"] += 1
            status.lifton_aa = v
            return ("", v)

        monkeypatch.setattr(lifton_class.Lifton_GENE, "orf_search_protein", _stub_orf)

        status = lifton_class.Lifton_Status()
        run_liftoff.process_liftoff_with_protein(
            SimpleNamespace(seqid="chr1"), gene, trans, {}, None, {}, None,
            "ref1", {}, {}, None, False, status, False)

        assert status.annotation == "LiftOn_miniprot"
        assert calls["n"] == 3
        # EXACTLY miniprot's 2 CDS-bearing exons -- no stale Liftoff-only exon,
        # no duplicates (the corruption class the chicken A/B caught).
        assert len(trans.exons) == 2
        coords = sorted((e.entry.start, e.entry.end) for e in trans.exons)
        assert coords == [(150, 250), (350, 450)]
        assert all(e.cds is not None for e in trans.exons)


# ═══════════════════════════════════════════════════════════════════════════
# The --no-miniprot-candidate CLI opt-out + resolver (Phase 1 formalization).
# resolve_miniprot_candidate_args publishes the CLI choice to run_liftoff's
# module default; the LIFTON_MINIPROT_CANDIDATE env var still WINS at read time
# (so the 24-cell matrix + A/B harness, which force via env, are unaffected).
# ═══════════════════════════════════════════════════════════════════════════
class TestResolveMiniprotCandidateArgs:
    @pytest.fixture(autouse=True)
    def _restore_module_default(self):
        # resolve_ mutates the module global; snapshot + restore so a False set
        # here can't leak into the default-ON behaviour tests above.
        orig = run_liftoff._MINIPROT_CANDIDATE_DEFAULT
        yield
        run_liftoff._MINIPROT_CANDIDATE_DEFAULT = orig

    def test_default_on(self, monkeypatch):
        monkeypatch.delenv("LIFTON_MINIPROT_CANDIDATE", raising=False)
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_candidate=True)
        lifton_main.resolve_miniprot_candidate_args(args)
        assert args.miniprot_candidate is True
        assert run_liftoff._MINIPROT_CANDIDATE_DEFAULT is True
        assert run_liftoff._miniprot_candidate_enabled() is True

    def test_no_flag_sets_default_off(self, monkeypatch):
        monkeypatch.delenv("LIFTON_MINIPROT_CANDIDATE", raising=False)
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_candidate=False)  # --no-miniprot-candidate
        lifton_main.resolve_miniprot_candidate_args(args)
        assert args.miniprot_candidate is False
        assert run_liftoff._MINIPROT_CANDIDATE_DEFAULT is False
        assert run_liftoff._miniprot_candidate_enabled() is False

    def test_env_wins_over_cli_default(self, monkeypatch):
        from lifton import lifton as lifton_main
        # CLI default ON, but env forces OFF/ON at read time.
        args = SimpleNamespace(miniprot_candidate=True)
        lifton_main.resolve_miniprot_candidate_args(args)
        monkeypatch.setenv("LIFTON_MINIPROT_CANDIDATE", "0")
        assert run_liftoff._miniprot_candidate_enabled() is False
        monkeypatch.setenv("LIFTON_MINIPROT_CANDIDATE", "1")
        assert run_liftoff._miniprot_candidate_enabled() is True

    def test_missing_attr_defaults_on(self, monkeypatch):
        monkeypatch.delenv("LIFTON_MINIPROT_CANDIDATE", raising=False)
        from lifton import lifton as lifton_main
        args = SimpleNamespace()  # no .miniprot_candidate attribute
        lifton_main.resolve_miniprot_candidate_args(args)
        assert args.miniprot_candidate is True
        assert run_liftoff._miniprot_candidate_enabled() is True
