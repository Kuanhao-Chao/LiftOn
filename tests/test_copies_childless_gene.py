"""A Liftoff ``-copies`` extra copy must keep its transcript, exons and CDS.

Reported by a user: with ``-copies``, LiftOn emitted only the *gene* line for an
extra gene copy, where Liftoff emits gene + mRNA + exon + CDS.

Root cause. Liftoff suffixes every feature of an extra copy with
``_<extra_copy_number>`` (``liftoff/write_new_gff.py:edit_copy_ids``), so the
copy arrives as ``gene-X_1`` / ``rna-X_1``. ``lifton_utils.get_ref_ids_liftoff``
resolved the GENE against the gene-keyed ``ref_features_dict`` (``gene-X_1`` ->
``gene-X``, correct), but for the TRANSCRIPT it called
``get_ID_base(liftoff_trans_id, None)`` -- and with ``None`` for the dict that
helper is the identity function. ``rna-X_1`` was therefore looked up verbatim,
``ref_db[...]`` raised, and ``run_liftoff.process_liftoff`` returned ``None``,
dropping the transcript and every exon and CDS under it while the already-built
gene stayed in the output as a bare line.

The conservative guard in ``get_ID_base`` is right -- reference ids that
genuinely end in ``_<int>`` exist (``...mrna.FMUND_1``) -- but passing ``None``
made it unconditional rather than conditional. ``resolve_ref_trans_id`` restores
the condition by checking the reference annotation itself, exact id first.

Measured on the 17-genome benchmark corpus before the fix: ~4,400 childless copy
genes (rice 539 of 815 copy genes, human->zebrafish 1,178 of 1,881), identical in
v1.0.11 and v1.0.12, i.e. shipping in every release. Every one of the 1,908 cases
checked across rice, human->zebrafish and arabidopsis resolved by stripping the
suffix, and in every case the resolved base was a child of the reference gene the
copy's own gene id resolved to.
"""

from __future__ import annotations

import textwrap

import pytest

from lifton import lifton_utils


# ---------------------------------------------------------------------------
# resolve_ref_trans_id — the rule in isolation
# ---------------------------------------------------------------------------

class _Feature:
    def __init__(self, feature_id, attributes=None):
        self.id = feature_id
        self.attributes = attributes or {}


class _FakeRefDb:
    """Minimal ``ref_db[id]`` stand-in; raises KeyError like _RefDbProxy."""

    def __init__(self, features):
        self._features = features

    def __getitem__(self, feature_id):
        return self._features[feature_id]


class TestResolveRefTransId:
    def test_copy_suffix_is_stripped_when_the_base_is_in_the_reference(self):
        db = _FakeRefDb({
            "rna-X": _Feature("rna-X", {"Parent": ["gene-X"]}),
        })
        assert lifton_utils.resolve_ref_trans_id(db, "rna-X_1", "gene-X") == "rna-X"

    def test_exact_id_wins_so_a_real_underscore_int_id_survives(self):
        """``FMUND_1`` is a real reference id, not a copy of ``FMUND``.

        Both ids exist here. Exact-first must return the id untouched; stripping
        would silently retarget the transcript onto a different feature.
        """
        db = _FakeRefDb({
            "FMUND": _Feature("FMUND", {"Parent": ["gene-F"]}),
            "FMUND_1": _Feature("FMUND_1", {"Parent": ["gene-F"]}),
        })
        assert lifton_utils.resolve_ref_trans_id(db, "FMUND_1", "gene-F") == "FMUND_1"

    def test_copy_base_belonging_to_another_gene_is_rejected(self):
        """The base exists but hangs off a different gene -> not this copy."""
        db = _FakeRefDb({
            "rna-X": _Feature("rna-X", {"Parent": ["gene-OTHER"]}),
        })
        assert lifton_utils.resolve_ref_trans_id(db, "rna-X_1", "gene-X") is None

    def test_guard_is_skipped_when_the_gene_is_unknown(self):
        """The 3-level hierarchy case resolves from the transcript id alone."""
        db = _FakeRefDb({
            "rna-X": _Feature("rna-X", {"Parent": ["gene-X"]}),
        })
        assert lifton_utils.resolve_ref_trans_id(db, "rna-X_1", None) == "rna-X"

    def test_unresolvable_id_returns_none_so_the_caller_still_warns(self):
        db = _FakeRefDb({})
        assert lifton_utils.resolve_ref_trans_id(db, "rna-X_1", "gene-X") is None
        assert lifton_utils.resolve_ref_trans_id(db, "rna-X", "gene-X") is None

    def test_non_numeric_suffix_is_not_a_copy_suffix(self):
        db = _FakeRefDb({"rna-X": _Feature("rna-X", {"Parent": ["gene-X"]})})
        assert lifton_utils.resolve_ref_trans_id(db, "rna-X_alt", "gene-X") is None

    def test_none_candidate(self):
        assert lifton_utils.resolve_ref_trans_id(_FakeRefDb({}), None, None) is None


class TestCopySuffixBase:
    @pytest.mark.parametrize("feature_id,expected", [
        ("rna-X_1", "rna-X"),
        ("rna-XM_015766011.2_1", "rna-XM_015766011.2"),   # real RefSeq shape
        ("gene-X_42", "gene-X"),
        ("rna-X", None),
        ("rna-X_alt", None),
        ("_1", None),            # no base
        ("", None),
        (None, None),
    ])
    def test_candidate_base(self, feature_id, expected):
        assert lifton_utils.copy_suffix_base(feature_id) == expected


def test_get_ref_ids_liftoff_is_unchanged_for_every_other_caller():
    """``get_ID_base(id, None)`` was the identity function.

    The both-ids branch used to return ``get_ID_base(liftoff_trans_id, None)``
    and now returns ``liftoff_trans_id``. Those are the same value for every
    input, so callers that resolve the id themselves (``run_evaluation``, via
    ``_lookup_reference_feature``) are untouched by this change.
    """
    from lifton.coreutils import get_ID_base
    for candidate in ("rna-X_1", "rna-X", "X", "0", "_1", "a_b_12",
                      "FMUND_1", "rna-XM_015766011.2_1", "", "x_"):
        assert get_ID_base(candidate, None) == candidate


# ---------------------------------------------------------------------------
# End to end, through the real pipeline
# ---------------------------------------------------------------------------

def _wrap(seq: str) -> str:
    return "\n".join(textwrap.wrap(seq, 60)) + "\n"


ORF = "ATG" + "GCT" * 31 + "TAA"          # 99 nt: M + A*31 + stop


@pytest.fixture
def hermetic_pipeline(monkeypatch):
    """Same contract as ``test_integration_pipeline``'s fixture: the external
    runners must never be reached, and the miniprot binary need not exist."""
    from lifton import lifton_utils as _lu, run_liftoff, run_miniprot

    monkeypatch.setattr(_lu, "check_miniprot_installed", lambda: None)
    monkeypatch.setattr(run_miniprot, "check_miniprot_installed", lambda: True)

    def _fail(*args, **kwargs):
        raise AssertionError(
            "External liftoff/miniprot must NOT be invoked in tests"
        )

    monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)
    monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)


@pytest.fixture
def copies_workspace(tmp_path):
    """A Liftoff output carrying an extra gene copy, in Liftoff's own shape.

    The copy is suffixed ``_1`` on the gene, the transcript, the exon and the
    CDS, and carries ``extra_copy_number=1`` on every row -- exactly what
    ``liftoff/write_new_gff.py`` writes. The target genome holds the same ORF at
    both loci so both copies lift cleanly.
    """
    work = tmp_path / "work"
    work.mkdir()

    chrom = ["A"] * 900
    for offset in (100, 500):                 # 1-based 101-199 and 501-599
        for i, ch in enumerate(ORF):
            chrom[offset + i] = ch
    seq = "".join(chrom)

    ref_fa = work / "ref.fa"
    ref_fa.write_text(">chr1\n" + _wrap(seq))
    tgt_fa = work / "tgt.fa"
    tgt_fa.write_text(">chr1\n" + _wrap(seq))

    # The reference knows only the single, unsuffixed locus.
    ref_gff = work / "ref.gff3"
    ref_gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
    )

    # Liftoff found two copies. Copy 0 keeps the plain ids; copy 1 is suffixed.
    liftoff_gff = work / "liftoff.gff3"
    liftoff_gff.write_text(
        "##gff-version 3\n"
        "chr1\tLiftoff\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding;extra_copy_number=0\n"
        "chr1\tLiftoff\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1;extra_copy_number=0\n"
        "chr1\tLiftoff\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1;extra_copy_number=0\n"
        "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1;extra_copy_number=0\n"
        "chr1\tLiftoff\tgene\t501\t599\t.\t+\t.\tID=gene1_1;gene_biotype=protein_coding;extra_copy_number=1\n"
        "chr1\tLiftoff\tmRNA\t501\t599\t.\t+\t.\tID=tx1_1;Parent=gene1_1;extra_copy_number=1\n"
        "chr1\tLiftoff\texon\t501\t599\t.\t+\t.\tID=exon1_1;Parent=tx1_1;extra_copy_number=1\n"
        "chr1\tLiftoff\tCDS\t501\t599\t.\t+\t0\tID=cds1_1;Parent=tx1_1;extra_copy_number=1\n"
    )

    # One miniprot hit over the primary locus, as in ``integration_workspace``.
    # It overlaps a lifted gene, so Step 8 suppresses it and it does not add a
    # model -- the copy hierarchy under test comes from Liftoff alone.
    miniprot_gff = work / "miniprot.gff3"
    miniprot_gff.write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t199\t.\t+\t.\tID=MP1;Target=tx1 1 33\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
    )

    out_dir = work / "out"
    out_dir.mkdir()
    return {
        "ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": ref_gff,
        "liftoff": liftoff_gff, "miniprot": miniprot_gff,
        "out": out_dir, "work": work,
    }


def _run(workspace, out_name, extra_argv=()):
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / out_name
    argv = [
        str(workspace["tgt_fa"]), str(workspace["ref_fa"]),
        "-g", str(workspace["ref_gff"]),
        "-L", str(workspace["liftoff"]),
        "-M", str(workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq", "--force",
        *extra_argv,
    ]
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    return out_gff


def _parse(path):
    """Return {feature_id: (featuretype, parent)} plus children-by-parent."""
    rows, children = {}, {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        attrs = dict(
            field.split("=", 1) for field in cols[8].split(";") if "=" in field
        )
        feature_id = attrs.get("ID")
        parent = attrs.get("Parent")
        if feature_id:
            rows[feature_id] = (cols[2], parent)
        if parent:
            children.setdefault(parent, []).append((cols[2], feature_id))
    return rows, children


class TestCopyKeepsItsHierarchy:
    def test_extra_copy_emits_transcript_exon_and_cds(
            self, copies_workspace, hermetic_pipeline):
        out_gff = _run(copies_workspace, "lifton.gff3")
        rows, children = _parse(out_gff)

        copy_genes = [
            fid for fid, (ftype, _) in rows.items()
            if ftype == "gene" and fid.endswith("_1")
        ]
        assert copy_genes, f"no copy gene was emitted at all: {sorted(rows)}"

        for gene_id in copy_genes:
            kinds = {ftype for ftype, _ in children.get(gene_id, [])}
            assert "mRNA" in kinds, (
                f"{gene_id} was emitted as a bare gene line -- the copy lost its "
                f"transcript. children={children.get(gene_id)}"
            )
            transcript_ids = [
                fid for ftype, fid in children[gene_id] if ftype == "mRNA"
            ]
            for transcript_id in transcript_ids:
                sub = {ftype for ftype, _ in children.get(transcript_id, [])}
                assert "exon" in sub and "CDS" in sub, (
                    f"{transcript_id} under {gene_id} has no {sub and 'CDS' or 'children'}"
                )

    def test_transcript_parent_matches_the_emitted_gene_id(
            self, copies_workspace, hermetic_pipeline):
        """The copy suffix on the gene comes from LiftOn's own counter, the one
        on the transcript from Liftoff's. The emitted Parent must still match."""
        out_gff = _run(copies_workspace, "lifton.gff3")
        rows, _ = _parse(out_gff)
        gene_ids = {fid for fid, (ftype, _) in rows.items() if ftype == "gene"}
        for fid, (ftype, parent) in rows.items():
            if ftype == "mRNA":
                assert parent in gene_ids, f"{fid} is orphaned (Parent={parent})"

    def test_both_copies_are_emitted(self, copies_workspace, hermetic_pipeline):
        out_gff = _run(copies_workspace, "lifton.gff3")
        rows, _ = _parse(out_gff)
        genes = [fid for fid, (ftype, _) in rows.items() if ftype == "gene"]
        mrnas = [fid for fid, (ftype, _) in rows.items() if ftype == "mRNA"]
        assert len(genes) == 2, f"expected both copies, got {genes}"
        assert len(mrnas) == 2, f"expected a transcript per copy, got {mrnas}"

    def test_no_childless_gene_is_reported(
            self, copies_workspace, hermetic_pipeline):
        """The run manifest must report zero genes emitted without children.

        This is the counter that would have made the bug visible: before the fix
        it was the number of extra copies, and nothing printed it.
        """
        import json

        out_gff = _run(copies_workspace, "lifton.gff3")
        manifest = json.loads(
            (out_gff.parent / "lifton_output" / "run_manifest.json").read_text()
        )
        assert manifest["counts"]["genes_emitted_without_children"] == 0


def test_parallel_matches_serial_on_a_copies_annotation(
        copies_workspace, hermetic_pipeline):
    """Pins the parallel prefetch.

    Workers read the reference through ``locus_pipeline._RefDbProxy``, which
    raises ``KeyError`` for any id the parent thread did not pre-fetch. The copy
    base is a second candidate id, so ``_populate_ref_attrs_for_descent`` has to
    cache it too -- otherwise the fix works at ``-t 1`` and silently does not at
    ``-t 4``, which is the default once threads are on.
    """
    serial = _run(copies_workspace, "serial.gff3", ("-t", "1"))
    parallel = _run(copies_workspace, "parallel.gff3",
                    ("-t", "4", "--locus-pipeline"))
    assert serial.read_text() == parallel.read_text()
