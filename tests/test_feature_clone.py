"""`clone_feature` must be indistinguishable from `copy.deepcopy` for LiftOn.

`copy.deepcopy` was the SECOND-hottest entry in both Step-7 profiles -- 36.2 M
calls on drosophila, 34.6 M on dog->cat, with 32.0 s cumulative inside
`run_liftoff._snapshot_merge_state` alone. Measured on a CDS carrying a full
RefSeq attribute set, a purpose-built clone is ~6.6x cheaper (33.1 us -> 5.0 us):
43% of the deepcopy went into recursing through attribute strings and 24% into
re-copying the `dialect` dict, which is a single shared parser-metadata object
that nothing in LiftOn ever writes.

Sharing `dialect` is the load-bearing assumption, so it is asserted here in both
directions: the clone shares it (that is the win) AND no LiftOn module writes to
it (that is what makes sharing safe).
"""
import copy
import tempfile
import os

import gffutils
import pytest

from lifton import coreutils


GFF3 = """##gff-version 3
chr1\tRefSeq\tgene\t1\t900\t.\t+\t.\tID=gene-A;Dbxref=GeneID:1,GenBank:X;gbkey=Gene;gene=A;locus_tag=L1
chr1\tRefSeq\tmRNA\t1\t900\t.\t+\t.\tID=rna-A;Parent=gene-A;gbkey=mRNA;gene=A;product=thing
chr1\tRefSeq\texon\t1\t300\t.\t+\t.\tID=exon-A-1;Parent=rna-A;gbkey=mRNA;gene=A;product=thing
chr1\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cds-A;Parent=rna-A;Dbxref=GeneID:1,GenBank:NP_1;Name=NP_1;gbkey=CDS;gene=A;locus_tag=L1;product=thing;protein_id=NP_1
"""

SCALAR_FIELDS = (
    "seqid", "source", "featuretype", "start", "end", "score", "strand",
    "frame", "id", "bin", "file_order", "keep_order", "sort_attribute_values",
)


@pytest.fixture(scope="module")
def db():
    path = tempfile.mktemp(suffix=".gff3")
    with open(path, "w") as handle:
        handle.write(GFF3)
    try:
        yield gffutils.create_db(path, ":memory:",
                                 merge_strategy="create_unique",
                                 keep_order=True)
    finally:
        os.unlink(path)


@pytest.fixture
def rich_cds(db):
    return db["cds-A"]


# ── equivalence with deepcopy ───────────────────────────────────────────────

@pytest.mark.parametrize("fid", ["gene-A", "rna-A", "exon-A-1", "cds-A"])
def test_every_field_matches_deepcopy(db, fid):
    source = db[fid]
    deep = copy.deepcopy(source)
    clone = coreutils.clone_feature(source)

    assert sorted(vars(clone)) == sorted(vars(deep))
    for field in SCALAR_FIELDS:
        assert getattr(clone, field) == getattr(deep, field), field
    assert dict(clone.attributes) == dict(deep.attributes)
    assert list(clone.extra) == list(deep.extra)
    assert clone.dialect == deep.dialect
    assert type(clone) is type(deep)
    assert type(clone.attributes) is type(source.attributes)


def test_in_memory_feature_with_a_plain_dict(make_gffutils_feature):
    source = make_gffutils_feature(
        featuretype="CDS", start=5, end=50, feature_id="c1",
        attributes={"ID": ["c1"], "Parent": ["t1"], "product": ["x"]})
    clone = coreutils.clone_feature(source)
    assert dict(clone.attributes) == dict(source.attributes)
    assert clone.start == 5 and clone.end == 50


def test_rendered_line_is_identical(rich_cds):
    from lifton.io import gff3_writer
    assert (gff3_writer.format_feature(coreutils.clone_feature(rich_cds))
            == gff3_writer.format_feature(rich_cds))


# ── independence ────────────────────────────────────────────────────────────

def test_mutating_the_clone_does_not_touch_the_source(rich_cds):
    clone = coreutils.clone_feature(rich_cds)
    clone.start = 999
    clone.attributes["ID"] = ["mutated"]
    clone.attributes["product"][0] = "mutated"
    clone.attributes["brand_new"] = ["x"]
    clone.extra.append("junk")

    assert rich_cds.start == 1
    assert rich_cds.attributes["ID"] == ["cds-A"]
    assert rich_cds.attributes["product"] == ["thing"]
    assert "brand_new" not in rich_cds.attributes
    assert list(rich_cds.extra) == []


def test_mutating_the_source_does_not_touch_the_clone(rich_cds):
    clone = coreutils.clone_feature(rich_cds)
    original = list(rich_cds.attributes["product"])
    try:
        rich_cds.attributes["product"] = ["changed"]
        assert clone.attributes["product"] == original
    finally:
        rich_cds.attributes["product"] = original


def test_attribute_lists_are_not_shared(rich_cds):
    clone = coreutils.clone_feature(rich_cds)
    for key, value in rich_cds.attributes.items():
        if isinstance(value, list):
            assert clone.attributes[key] is not value, key


# ── the OTHER Feature class ─────────────────────────────────────────────────

class TestGffbaseSlotsFeature:
    """`lifton.gffbase.feature.Feature` declares __slots__ and has NO __dict__.

    It is the class the `--inmemory-liftoff` path uses, and a __dict__-only
    clone silently produced a different lift there -- caught by the 24-cell
    matrix, not by any unit test, which is exactly why this class is covered
    here now.
    """

    @staticmethod
    def _feature():
        from lifton.gffbase.feature import Feature
        return Feature(
            seqid="chr1", source="Liftoff", featuretype="CDS", start=10,
            end=99, score=".", strand="+", frame="0",
            attributes={"ID": ["cds-1"], "Parent": ["rna-1"],
                        "product": ["thing"]},
            id="cds-1",
        )

    def test_the_class_really_has_no_dict(self):
        feature = self._feature()
        assert getattr(feature, "__dict__", None) is None

    def test_clone_matches_deepcopy_field_for_field(self):
        feature = self._feature()
        deep = copy.deepcopy(feature)
        clone = coreutils.clone_feature(feature)
        for name in coreutils._iter_slots(type(feature)):
            if name in ("attributes", "extra", "_attributes_blob"):
                continue
            assert getattr(clone, name, None) == getattr(deep, name, None), name
        assert dict(clone.attributes) == dict(deep.attributes)

    def test_clone_is_independent(self):
        feature = self._feature()
        clone = coreutils.clone_feature(feature)
        clone.source = "LiftOn"
        clone.attributes["product"] = ["mutated"]
        assert feature.source == "Liftoff"
        assert feature.attributes["product"] == ["thing"]

    def test_unset_slots_do_not_raise(self):
        from lifton.gffbase.feature import Feature
        bare = Feature.__new__(Feature)
        bare.seqid = "chr1"
        bare.attributes = {}
        clone = coreutils.clone_feature(bare)
        assert clone.seqid == "chr1"


# ── the dialect decision, asserted in both directions ───────────────────────

def test_dialect_is_shared_on_purpose(rich_cds):
    clone = coreutils.clone_feature(rich_cds)
    assert clone.dialect is rich_cds.dialect


def test_a_database_already_shares_one_dialect(db):
    features = list(db.all_features())
    assert len(features) > 1
    assert len({id(f.dialect) for f in features}) == 1


def test_no_lifton_module_writes_to_a_dialect():
    """Sharing is only safe while `dialect` stays read-only.

    A write would look like `feature.dialect[...] = ...` or a mutating method
    call on it. Vendored Liftoff only reads `dialect['fmt']`.
    """
    import pathlib
    import re

    root = pathlib.Path(coreutils.__file__).parent
    writes = []
    pattern = re.compile(
        r"\.dialect\s*\[[^\]]*\]\s*=|"
        r"\.dialect\.(?:update|pop|clear|setdefault|__setitem__)\s*\("
    )
    for path in root.rglob("*.py"):
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            if pattern.search(line):
                writes.append(f"{path.relative_to(root)}:{lineno}")
    assert writes == [], f"dialect is written here, sharing is unsafe: {writes}"


# ── clone_attributes on its own ─────────────────────────────────────────────

def test_clone_attributes_matches_deepcopy(rich_cds):
    source = rich_cds.attributes
    assert dict(coreutils.clone_attributes(source)) == dict(copy.deepcopy(source))


def test_clone_attributes_passes_none_through():
    assert coreutils.clone_attributes(None) is None


def test_clone_attributes_keeps_scalar_values():
    cloned = coreutils.clone_attributes({"a": "plain", "b": ["l"], "c": ("t",)})
    assert cloned["a"] == "plain"
    assert cloned["b"] == ["l"] and isinstance(cloned["b"], list)
    assert cloned["c"] == ["t"]


# ── the Lifton_* __deepcopy__ overrides ─────────────────────────────────────

class TestLiftonDeepcopy:
    @staticmethod
    def _exon(db):
        from lifton.lifton_class import Lifton_EXON
        exon = Lifton_EXON(coreutils.clone_feature(db["exon-A-1"]))
        exon.add_cds(coreutils.clone_feature(db["cds-A"]))
        return exon

    def test_exon_clone_is_independent(self, db):
        exon = self._exon(db)
        clone = copy.deepcopy(exon)

        assert clone is not exon and clone.cds is not exon.cds
        assert clone.entry.start == exon.entry.start
        assert dict(clone.cds.entry.attributes) == dict(exon.cds.entry.attributes)

        clone.entry.start = 12345
        clone.cds.entry.attributes["product"] = ["mutated"]
        assert exon.entry.start != 12345
        assert exon.cds.entry.attributes["product"] != ["mutated"]

    def test_exon_without_a_cds(self, db):
        from lifton.lifton_class import Lifton_EXON
        exon = Lifton_EXON(coreutils.clone_feature(db["exon-A-1"]))
        assert copy.deepcopy(exon).cds is None

    def test_cds_clone_preserves_source_id_provenance(self, db):
        from lifton.lifton_class import Lifton_CDS
        entry = coreutils.clone_feature(db["cds-A"])
        entry.attributes["extra_copy_number"] = ["2"]
        entry.attributes["ID"] = ["cds-A_2"]
        cds = Lifton_CDS(entry)
        assert cds._source_id == "cds-A_2" and cds._source_id_base == "cds-A"

        clone = copy.deepcopy(cds)
        # Re-running __init__ on the clone would lose this: __init__ derives the
        # base from `extra_copy_number`, which it then POPS.
        assert clone._source_id == "cds-A_2"
        assert clone._source_id_base == "cds-A"
        assert clone._source_copy_number == "2"

    def test_memo_keeps_a_shared_exon_shared(self, db):
        exon = self._exon(db)
        pair = copy.deepcopy([exon, exon])
        assert pair[0] is pair[1]

    def test_deepcopying_the_exon_list_matches_the_old_behaviour(self, db):
        """`_snapshot_merge_state` deep-copies `trans.exons`; the snapshot must
        still restore coordinates and attributes exactly."""
        exons = [self._exon(db) for _ in range(3)]
        for offset, exon in enumerate(exons):
            exon.entry.start = 100 + offset
        snapshot = copy.deepcopy(exons)

        # The trial merge then rewrites the live exons, as update_cds_list does.
        for exon in exons:
            exon.entry.start = 999
            exon.cds = None

        assert [e.entry.start for e in snapshot] == [100, 101, 102]
        assert all(e.cds is not None for e in snapshot)
        assert all(e.cds.entry.attributes["product"] == ["thing"]
                   for e in snapshot)
