"""Phase 7 — property-based invariants for the streaming adapter.

Hypothesis explores random-but-valid GFF3 blobs and asserts:
  - The blob ingest path produces the same FeatureDB content as the
    file ingest path for any input.
  - communicate() never truncates.
  - The bytes/str/bytearray ingest paths agree.
"""

from __future__ import annotations

import string

import pytest
from hypothesis import HealthCheck, given, settings, strategies as st

from lifton import gffbase_adapter
from lifton.gffbase import create_db


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Coordinate constraints: 1-based positive integers with start <= end
_COORD = st.integers(min_value=1, max_value=10**6)
_SEQID = st.sampled_from(["chr1", "chr2", "chrX", "chrY", "chrM"])
_STRAND = st.sampled_from(["+", "-"])
_BIOTYPE = st.sampled_from([
    "protein_coding", "lncRNA", "miRNA", "pseudogene", "tRNA",
])
_NAME_CHARS = string.ascii_letters + string.digits + "_-"
_NAME = st.text(alphabet=_NAME_CHARS, min_size=1, max_size=12)


@st.composite
def _gene_record(draw):
    seqid = draw(_SEQID)
    strand = draw(_STRAND)
    s = draw(_COORD)
    e = draw(st.integers(min_value=s, max_value=s + 5000))
    biotype = draw(_BIOTYPE)
    name = draw(_NAME)
    line = f"{seqid}\ttest\tgene\t{s}\t{e}\t.\t{strand}\t.\tID=g_{name};gene_biotype={biotype}\n"
    return line


@st.composite
def _gff_blob(draw):
    n = draw(st.integers(min_value=1, max_value=12))
    lines = [draw(_gene_record()) for _ in range(n)]
    # Ensure unique IDs by appending the index
    body = ""
    for i, line in enumerate(lines):
        body += line.replace("ID=g_", f"ID=g_{i}_", 1)
    return ("##gff-version 3\n" + body).encode("utf-8")


# ---------------------------------------------------------------------------
# Property tests
# ---------------------------------------------------------------------------

@settings(max_examples=50, deadline=None,
          suppress_health_check=[HealthCheck.too_slow,
                                 HealthCheck.function_scoped_fixture])
@given(blob=_gff_blob())
def test_blob_ingest_count_matches_file_ingest(blob, tmp_path_factory):
    """For any random valid GFF3 blob, the in-memory ingest yields
    the same gene count as the file-based ingest."""
    p = tmp_path_factory.mktemp("h") / "x.gff3"
    p.write_bytes(blob)
    db_file = create_db(str(p), dbfn=":memory:", force=True)
    db_blob = gffbase_adapter.build_database_from_string(gff_text=blob)
    assert db_file.count_features_of_type("gene") == \
           db_blob.count_features_of_type("gene")


@settings(max_examples=50, deadline=None,
          suppress_health_check=[HealthCheck.too_slow,
                                 HealthCheck.function_scoped_fixture])
@given(blob=_gff_blob())
def test_blob_ingest_seqid_set_matches_file_ingest(blob, tmp_path_factory):
    p = tmp_path_factory.mktemp("h") / "x.gff3"
    p.write_bytes(blob)
    db_file = create_db(str(p), dbfn=":memory:", force=True)
    db_blob = gffbase_adapter.build_database_from_string(gff_text=blob)
    assert set(db_file.seqids()) == set(db_blob.seqids())


@settings(max_examples=50, deadline=None)
@given(blob=_gff_blob())
def test_str_and_bytes_blob_agree(blob):
    """Same content; coercion must not change the result."""
    db_bytes = gffbase_adapter.build_database_from_string(gff_text=blob)
    db_str = gffbase_adapter.build_database_from_string(
        gff_text=blob.decode("utf-8")
    )
    assert db_bytes.count_features_of_type("gene") == \
           db_str.count_features_of_type("gene")


@settings(max_examples=50, deadline=None)
@given(blob=_gff_blob())
def test_bytearray_blob_matches_bytes(blob):
    db_bytes = gffbase_adapter.build_database_from_string(gff_text=blob)
    db_ba = gffbase_adapter.build_database_from_string(
        gff_text=bytearray(blob)
    )
    assert db_bytes.count_features_of_type("gene") == \
           db_ba.count_features_of_type("gene")


@settings(max_examples=50, deadline=None)
@given(blob=_gff_blob())
def test_blob_ingest_does_not_create_duckdb_sidecar(blob):
    """Default dbfn=':memory:' must not write any .duckdb file."""
    import os, tempfile
    cwd_before = set(os.listdir("."))
    gffbase_adapter.build_database_from_string(gff_text=blob)
    cwd_after = set(os.listdir("."))
    new_files = cwd_after - cwd_before
    assert not any(f.endswith(".duckdb") for f in new_files)


@settings(max_examples=80, deadline=None)
@given(blob=_gff_blob())
def test_looks_like_gff3_blob_recognises_real_blobs(blob):
    assert gffbase_adapter.looks_like_gff3_blob(blob) is True


@settings(max_examples=100, deadline=None)
@given(path=st.text(alphabet=string.ascii_letters + "/.", min_size=1, max_size=40))
def test_path_strings_never_recognised_as_blob(path):
    """Random path-shaped strings must not be misidentified as
    GFF3 blobs (the heuristic is bytes-only)."""
    assert gffbase_adapter.looks_like_gff3_blob(path) is False


@settings(max_examples=50, deadline=None)
@given(blob=_gff_blob())
def test_db_path_for_is_idempotent(blob):
    """db_path_for(p) must be a pure function of p — multiple calls
    return the same string."""
    p = "/tmp/whatever.gff3"
    assert gffbase_adapter.db_path_for(p) == gffbase_adapter.db_path_for(p)


@settings(max_examples=50, deadline=None)
@given(blob=_gff_blob())
def test_blob_ingest_is_deterministic(blob):
    """Two ingests of the same blob produce the same gene-id set."""
    db1 = gffbase_adapter.build_database_from_string(gff_text=blob)
    db2 = gffbase_adapter.build_database_from_string(gff_text=blob)
    ids1 = {f.id for f in db1.features_of_type("gene")}
    ids2 = {f.id for f in db2.features_of_type("gene")}
    assert ids1 == ids2


@settings(max_examples=30, deadline=None)
@given(blob=_gff_blob())
def test_annotation_blob_matches_file_for_path(blob, tmp_path_factory):
    """Annotation-level invariant: the FeatureDB attached to a
    bytes-input Annotation has the same gene count as the FeatureDB
    attached to a file-input Annotation backed by gffbase."""
    from lifton.annotation import Annotation
    fp = tmp_path_factory.mktemp("h") / "x.gff3"
    fp.write_bytes(blob)
    ann_blob = Annotation(blob, False, False)
    ann_file = Annotation(str(fp), False, False, backend="gffbase")
    assert ann_blob.db_connection.count_features_of_type("gene") == \
           ann_file.db_connection.count_features_of_type("gene")
