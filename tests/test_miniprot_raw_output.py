"""Focused tests for miniprot raw-output capture and bounded aggregation."""

from io import BytesIO
from types import SimpleNamespace
from unittest import mock

from lifton import run_miniprot
from lifton.native_bindings import GFF3Hit, MiniprotIndex
from lifton.native_bindings import miniprot_facade


GFF_BLOB = (
    b"##gff-version 3\n"
    b"chr1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
    b"ID=MP1;Target=tx1 1 30\n"
    b"chr1\tminiprot\tCDS\t1\t90\t.\t+\t0\t"
    b"ID=MP1.cds1;Parent=MP1\n"
)


def _facade_process(stdout=GFF_BLOB, stderr=b"", returncode=0):
    proc = mock.MagicMock()
    proc.communicate.return_value = (stdout, stderr)
    proc.returncode = returncode
    proc.wait.return_value = returncode
    proc.stdout = BytesIO(stdout)
    proc.stderr = BytesIO(stderr)
    return proc


def test_align_all_raw_only_skips_parser_and_preserves_bytes():
    """Raw mode must not decode/materialize a duplicate hit graph."""
    idx = MiniprotIndex("target.fa", ref_proteins_path="proteins.fa")

    with mock.patch.object(
        miniprot_facade.subprocess,
        "Popen",
        return_value=_facade_process(),
    ), mock.patch.object(
        GFF3Hit,
        "from_gff_line",
        side_effect=AssertionError("raw-only mode invoked the parser"),
    ) as parser:
        bundle = idx.align_all(raw_only=True)

    assert bundle.hits == []
    assert bundle.raw_bytes == GFF_BLOB
    assert idx.raw_bytes == GFF_BLOB
    parser.assert_not_called()


def test_parsed_call_reuses_bytes_cached_by_raw_only_call():
    """Opting into raw mode must not break the facade's default contract."""
    idx = MiniprotIndex("target.fa", ref_proteins_path="proteins.fa")
    proc = _facade_process()

    with mock.patch.object(
        miniprot_facade.subprocess,
        "Popen",
        return_value=proc,
    ) as popen:
        raw_bundle = idx.align_all(raw_only=True)
        parsed_bundle = idx.align_all()

    assert raw_bundle.raw_bytes == parsed_bundle.raw_bytes == GFF_BLOB
    assert [hit.featuretype for hit in parsed_bundle.hits] == ["mRNA", "CDS"]
    popen.assert_called_once()


def test_real_binding_shape_raw_only_skips_reparse_and_preserves_lines():
    """The forward-compatible binding branch receives the same fast path."""
    lines = [
        "chr2\tminiprot\tmRNA\t5\t100\t.\t-\t.\tID=MP2",
        "chr2\tminiprot\tCDS\t5\t100\t.\t-\t0\tParent=MP2",
    ]
    native_hits = [
        SimpleNamespace(to_gff_line=lambda line=line: line)
        for line in lines
    ]
    native = SimpleNamespace(align_file=lambda proteins: iter(native_hits))
    idx = MiniprotIndex("target.fa", ref_proteins_path="proteins.fa")
    idx._native = native

    with mock.patch.object(
        GFF3Hit,
        "from_gff_line",
        side_effect=AssertionError("native raw-only mode invoked the parser"),
    ) as parser:
        bundle = idx.align_all(raw_only=True)

    assert bundle.raw_bytes == ("\n".join(lines) + "\n").encode("utf-8")
    assert bundle.hits == []
    parser.assert_not_called()


def test_native_stream_pipeline_writes_database_without_hit_graph(tmp_path):
    """Native flag keeps the bounded DB stream path parser-free."""
    args = SimpleNamespace(
        mp_options="",
        stream=True,
        native=True,
        threads=1,
    )

    with mock.patch.object(
        miniprot_facade.subprocess,
        "Popen",
        return_value=_facade_process(),
    ), mock.patch.object(
        GFF3Hit,
        "from_gff_line",
        side_effect=AssertionError("pipeline invoked the parser"),
    ) as parser:
        output = run_miniprot.run_miniprot(
            str(tmp_path) + "/", args, "target.fa", "proteins.fa",
        )

    assert output.endswith("miniprot.duckdb")
    from lifton.gffbase_adapter import open_database_path
    db = open_database_path(output)
    assert db.count_features_of_type("mRNA") == 1
    db.conn.close()
    parser.assert_not_called()


def test_chunk_drain_preserves_streams_across_many_reads():
    """BytesIO aggregation is byte-exact for both concurrently drained pipes."""
    stdout = b"abcdefghijklmno"
    stderr = b"warning: one\nwarning: two\n"
    proc = SimpleNamespace(
        stdout=BytesIO(stdout),
        stderr=BytesIO(stderr),
        wait=lambda: 7,
    )

    captured_out, captured_err, returncode = run_miniprot._drain_stream_chunks(
        proc, chunk_size=3,
    )

    assert captured_out == stdout
    assert captured_err == stderr
    assert returncode == 7
