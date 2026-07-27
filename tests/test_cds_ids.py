"""GitHub issue #32 / #8 — CDS features must carry an ``ID=``.

Reconstructed and novel CDS paths historically reset each CDS's attributes to
``{Parent}`` only, so emitted CDS lines had no ``ID=``.
``Lifton_TRANS.normalize_containment`` (the write funnel, Iter-24) preserves a
stable source CDS ID when present and otherwise assigns a shared
``cds-<trans>`` to a transcript's CDS segments — the valid GFF3
discontinuous-CDS form (same ID + same Parent + same type), which
``gff3_validator.DISCONTINUOUS_FEATURE_TYPES = {"CDS"}`` already anticipates.
Gated by the same ``LIFTON_NO_CONTAINMENT_NORMALIZE`` escape hatch as the rest of
the write-funnel normalization, so it reproduces the pre-fix bytes when disabled.
"""
import gffutils

from lifton.cds_id_allocator import CdsIdAllocator
from lifton.gff3_validator import validate_gff3_structure
from lifton.io.gff3_writer import format_feature
from lifton.lifton import _allocator_source_ids
from lifton.lifton_class import Lifton_TRANS, Lifton_EXON, Lifton_CDS


def _ent(ftype, s, e, fid, parent):
    f = gffutils.Feature(seqid="chr1", source="LiftOn", featuretype=ftype, start=s,
                         end=e, strand="+", frame=".",
                         attributes={"ID": [fid], "Parent": [parent]})
    f.id = fid                       # the gffutils DB populates .id; set it here
    return f


def _trans(trans_id, cds_spans):
    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = _ent("mRNA", 101, 399, trans_id, "gene9")
    trans.exons = []
    for i, (cs, ce) in enumerate(cds_spans, 1):
        ex = Lifton_EXON(_ent("exon", cs, ce, f"exon-A-{i}", trans_id))
        ex.cds = Lifton_CDS(_ent("CDS", cs, ce, "stale", trans_id))
        ex.cds._source_id = None
        ex.cds._source_id_base = None
        ex.cds.entry.attributes = {"Parent": [trans_id]}   # add_lifton_cds resets (no ID)
        trans.exons.append(ex)
    return trans


def _set_source_id(cds, identifier):
    cds._source_id = identifier
    cds._source_id_base = identifier
    cds.entry.attributes["ID"] = [identifier]


def _write_transcripts(path, transcripts):
    gene = _ent("gene", 1, 500, "gene9", "")
    gene.attributes.pop("Parent", None)
    rows = ["##gff-version 3", format_feature(gene)]
    for trans in transcripts:
        rows.append(format_feature(trans.entry))
        rows.extend(format_feature(exon.entry) for exon in trans.exons)
        for exon in trans.exons:
            if exon.cds is not None:
                exon.cds.entry.frame = "0"
                rows.append(format_feature(exon.cds.entry))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_cds_segments_share_one_cds_id():
    trans = _trans("rna-NM_9", [(101, 199), (301, 399)])
    trans.normalize_containment()
    ids = [e.cds.entry.attributes["ID"] for e in trans.exons]
    assert ids == [["cds-NM_9"], ["cds-NM_9"]]        # shared across segments, "rna-" stripped


def test_cds_id_without_rna_prefix():
    trans = _trans("mytx", [(101, 199)])
    trans.normalize_containment()
    assert trans.exons[0].cds.entry.attributes["ID"] == ["cds-mytx"]


def test_existing_cds_id_is_preserved_and_fills_missing_segments():
    trans = _trans("rna-NM_9", [(101, 199), (301, 399)])
    _set_source_id(trans.exons[0].cds, "cds-NP_9.1")

    trans.normalize_containment()

    ids = [e.cds.entry.attributes["ID"] for e in trans.exons]
    assert ids == [["cds-NP_9.1"], ["cds-NP_9.1"]]


def test_reconstructed_cds_retains_source_id_provenance():
    exon = Lifton_EXON(_ent("exon", 101, 199, "exon-A-1", "rna-NM_9"))
    cds = Lifton_CDS(_ent("CDS", 101, 199, "cds-NP_9.1", "old-parent"))

    exon.add_lifton_cds(cds)

    assert exon.cds.entry.attributes == {"Parent": ["rna-NM_9"]}
    assert exon.cds._source_id == "cds-NP_9.1"


def test_inconsistent_existing_cds_ids_are_not_collapsed():
    trans = _trans("rna-NM_9", [(101, 149), (201, 249), (301, 399)])
    _set_source_id(trans.exons[0].cds, "cds-segment-1")
    _set_source_id(trans.exons[1].cds, "cds-segment-2")

    trans.normalize_containment()

    ids = [e.cds.entry.attributes["ID"] for e in trans.exons]
    assert ids == [["cds-segment-1"], ["cds-segment-2"], ["cds-NM_9"]]


def test_existing_cds_id_gets_copy_suffix():
    trans = _trans("rna-NM_9_1", [(101, 199), (301, 399)])
    trans.ref_tran_id = "rna-NM_9"
    for exon in trans.exons:
        _set_source_id(exon.cds, "cds-NP_9.1")

    trans.normalize_containment()

    ids = [e.cds.entry.attributes["ID"] for e in trans.exons]
    assert ids == [["cds-NP_9.1_1"], ["cds-NP_9.1_1"]]


def test_natural_numeric_suffix_does_not_mask_copy_suffix():
    trans = _trans("rna-source_1_1", [(101, 199)])
    trans.ref_tran_id = "rna-source_1"
    _set_source_id(trans.exons[0].cds, "stable_1")

    trans.normalize_containment()

    assert trans.exons[0].cds.entry.attributes["ID"] == ["stable_1_1"]


def test_liftoff_copy_suffix_is_not_appended_twice():
    trans = _trans("rna-NM_9_1", [(101, 199)])
    trans.ref_tran_id = "rna-NM_9"
    _set_source_id(trans.exons[0].cds, "cds-NP_9.1_1")
    trans.exons[0].cds._source_copy_number = "1"
    trans.exons[0].cds._source_id_base = "cds-NP_9.1"

    trans.normalize_containment()

    assert trans.exons[0].cds.entry.attributes["ID"] == ["cds-NP_9.1_1"]


def test_cds_id_colliding_with_exon_id_is_replaced():
    trans = _trans("rna-NM_9", [(101, 199)])
    trans.exons[0].entry.id = "MP9.cds1"
    trans.exons[0].entry.attributes["ID"] = ["MP9.cds1"]
    _set_source_id(trans.exons[0].cds, "MP9.cds1")

    trans.normalize_containment()

    assert trans.exons[0].entry.attributes["ID"] == ["MP9.cds1"]
    assert trans.exons[0].cds.entry.attributes["ID"] == ["cds-NM_9"]


def test_synthetic_cds_id_colliding_with_exon_id_is_disambiguated():
    trans = _trans("rna-NM_9", [(101, 199)])
    trans.exons[0].entry.id = "cds-NM_9"
    trans.exons[0].entry.attributes["ID"] = ["cds-NM_9"]

    trans.normalize_containment()

    assert trans.exons[0].cds.entry.attributes["ID"] == ["cds-NM_9-1"]


def test_lifton_cds_remembers_removed_source_copy_number():
    entry = _ent("CDS", 101, 199, "cds-NP_9.1_1", "rna-NM_9_1")
    entry.attributes["extra_copy_number"] = ["1"]

    cds = Lifton_CDS(entry)

    assert cds._source_id == "cds-NP_9.1_1"
    assert cds._source_copy_number == "1"
    assert "extra_copy_number" not in cds.entry.attributes


def test_file_allocator_reserves_future_stable_ids(tmp_path):
    allocator = CdsIdAllocator({"cds-X"})
    synthetic = _trans("rna-X", [(101, 199)])
    stable = _trans("rna-Y", [(301, 399)])
    stable.exons[0].entry.id = "exon-B-1"
    stable.exons[0].entry.attributes["ID"] = ["exon-B-1"]
    _set_source_id(stable.exons[0].cds, "cds-X")

    synthetic.normalize_containment(cds_id_allocator=allocator)
    stable.normalize_containment(cds_id_allocator=allocator)

    assert synthetic.exons[0].cds.entry.attributes["ID"] == ["cds-X-1"]
    assert stable.exons[0].cds.entry.attributes["ID"] == ["cds-X"]
    output = tmp_path / "collision-free.gff3"
    _write_transcripts(output, [synthetic, stable])
    assert validate_gff3_structure(str(output)).is_valid


def test_copy_cannot_claim_another_sources_exact_stable_id(tmp_path):
    allocator = CdsIdAllocator({"cds-A", "cds-A_1"})
    copied = _trans("rna-copy_1", [(101, 199)])
    copied.ref_tran_id = "rna-copy"
    copied_cds = copied.exons[0].cds
    _set_source_id(copied_cds, "cds-A_1")
    copied_cds._source_copy_number = "1"
    copied_cds._source_id_base = "cds-A"

    future_stable = _trans("rna-future", [(301, 399)])
    future_stable.exons[0].entry.id = "exon-future-1"
    future_stable.exons[0].entry.attributes["ID"] = ["exon-future-1"]
    _set_source_id(future_stable.exons[0].cds, "cds-A_1")

    copied.normalize_containment(cds_id_allocator=allocator)
    future_stable.normalize_containment(cds_id_allocator=allocator)

    copied_id = copied_cds.entry.attributes["ID"][0]
    future_id = future_stable.exons[0].cds.entry.attributes["ID"][0]
    assert copied_id != "cds-A_1"
    assert future_id == "cds-A_1"
    output = tmp_path / "future-stable-id.gff3"
    _write_transcripts(output, [copied, future_stable])
    assert validate_gff3_structure(str(output)).is_valid


def test_copy_checks_non_cds_source_ids_without_materializing_namespace():
    source_ids = {"stable", "stable_1"}
    allocator = CdsIdAllocator(
        exact_reserved_source_ids=source_ids,
    )
    copied = _trans("rna-copy_1", [(101, 199)])
    copied.ref_tran_id = "rna-copy"
    copied_cds = copied.exons[0].cds
    _set_source_id(copied_cds, "stable_1")
    copied_cds._source_copy_number = "1"
    copied_cds._source_id_base = "stable"
    future_stable = _trans("rna-future", [(301, 399)])
    _set_source_id(future_stable.exons[0].cds, "stable_1")

    copied.normalize_containment(cds_id_allocator=allocator)
    future_stable.normalize_containment(cds_id_allocator=allocator)

    assert copied_cds.entry.attributes["ID"] != ["stable_1"]
    assert future_stable.exons[0].cds.entry.attributes["ID"] == ["stable_1"]


def test_database_fallback_uses_declared_ids_not_internal_row_keys(tmp_path):
    source = tmp_path / "row-keys.gff3"
    database_path = tmp_path / "row-keys.db"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=gene1\n"
        "chr1\tRefSeq\tmRNA\t1\t200\t.\t+\t.\tID=rna1;Parent=gene1\n"
        "chr1\tRefSeq\tCDS\t1\t50\t.\t+\t0\tID=stable;Parent=rna1\n"
        "chr1\tRefSeq\tCDS\t101\t150\t.\t+\t0\tID=stable;Parent=rna1\n"
        "chr1\tRefSeq\tmRNA\t301\t400\t.\t+\t.\tID=rna2;Parent=gene1\n"
        "chr1\tRefSeq\tCDS\t301\t350\t.\t+\t0\t"
        "ID=stable_1;Parent=rna2\n",
        encoding="utf-8",
    )
    gffutils.create_db(
        str(source),
        dbfn=str(database_path),
        force=True,
        merge_strategy="create_unique",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
    database = gffutils.FeatureDB(str(database_path))

    assert database["stable_1"].attributes["ID"] == ["stable"]
    cds_namespace_ids, copy_suffix_ids = _allocator_source_ids(database)

    assert cds_namespace_ids == set()
    assert copy_suffix_ids == {"stable_1"}


def test_original_and_copy_with_natural_suffix_validate(tmp_path):
    allocator = CdsIdAllocator()
    original = _trans("rna-source_1", [(101, 199)])
    original.ref_tran_id = "rna-source_1"
    _set_source_id(original.exons[0].cds, "stable_1")
    copied = _trans("rna-source_1_1", [(301, 399)])
    copied.ref_tran_id = "rna-source_1"
    copied.exons[0].entry.id = "exon-copy-1"
    copied.exons[0].entry.attributes["ID"] = ["exon-copy-1"]
    _set_source_id(copied.exons[0].cds, "stable_1_1")
    copied.exons[0].cds._source_copy_number = "1"
    copied.exons[0].cds._source_id_base = "stable_1"

    original.normalize_containment(cds_id_allocator=allocator)
    copied.normalize_containment(cds_id_allocator=allocator)

    assert original.exons[0].cds.entry.attributes["ID"] == ["stable_1"]
    assert copied.exons[0].cds.entry.attributes["ID"] == ["stable_1_1"]
    output = tmp_path / "copy-safe.gff3"
    _write_transcripts(output, [original, copied])
    assert validate_gff3_structure(str(output)).is_valid


def test_normalization_is_idempotent_for_synthetic_and_preserved_ids():
    allocator = CdsIdAllocator()
    synthetic = _trans("rna-X", [(101, 199)])
    stable = _trans("rna-Y", [(301, 399)])
    _set_source_id(stable.exons[0].cds, "stable")

    for trans in (synthetic, stable):
        trans.normalize_containment(cds_id_allocator=allocator)
        first = trans.exons[0].cds.entry.attributes["ID"]
        trans.normalize_containment(cds_id_allocator=allocator)
        assert trans.exons[0].cds.entry.attributes["ID"] == first


def test_cds_id_reverts_under_escape_hatch(monkeypatch):
    monkeypatch.setenv("LIFTON_NO_CONTAINMENT_NORMALIZE", "1")
    trans = _trans("rna-NM_9", [(101, 199), (301, 399)])
    trans.normalize_containment()
    # normalize is a no-op -> CDS keep their ID-less {Parent} attributes.
    assert "ID" not in trans.exons[0].cds.entry.attributes


def test_reconstructed_source_id_is_not_materialized_under_escape_hatch(
        monkeypatch):
    monkeypatch.setenv("LIFTON_NO_CONTAINMENT_NORMALIZE", "1")
    trans = _trans("rna-NM_9", [(101, 199)])
    source = Lifton_CDS(
        _ent("CDS", 101, 199, "cds-NP_9.1", "old-parent")
    )
    trans.exons[0].add_lifton_cds(source)

    trans.normalize_containment()

    assert trans.exons[0].cds._source_id == "cds-NP_9.1"
    assert trans.exons[0].cds.entry.attributes == {
        "Parent": ["rna-NM_9"]
    }
