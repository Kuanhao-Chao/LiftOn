"""Focused regressions for lazy Step 7 and direct miniprot ingestion."""

from __future__ import annotations

import io
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import pytest

from lifton import locus_pipeline, parallel, run_miniprot
from lifton.gffbase import GFF3ChunkDecoder, create_db
from lifton.gffbase import ingest as gffbase_ingest
from lifton.gffbase import interface as gffbase_interface
from lifton.locus_pipeline import LocusResult, MaterialisedLocus, StepContext
from lifton.native_bindings import MiniprotIndex
from lifton.native_bindings import miniprot_facade


GFF = (
    b"##gff-version 3\n"
    b"chr1\tmp\tmRNA\t1\t90\t.\t+\t.\tID=MP1;Target=tx1 1 30\n"
    b"chr1\tmp\tCDS\t1\t90\t.\t+\t0\tID=CDS1;Parent=MP1\n"
)


class _Gene:
    ref_gene_id = "ref"

    def __init__(self, index):
        self.index = index

    def write_entry(self, handle, _stats):
        handle.write(f"{self.index}\n")
        return True


def _ctx(db, window=4):
    return StepContext(
        ref_db={}, l_feature_db=db, m_feature_db=None,
        ref_id_2_m_id_trans_dict={}, tree_dict={}, tgt_fai=None,
        ref_proteins={}, ref_trans={}, ref_features_dict={},
        fw_score=None, fw_chain=None,
        args=SimpleNamespace(
            native=False, writer_max_pending=0,
            step7_max_inflight=window, threads=3,
        ),
    )


def test_step7_root_scan_is_lazy_and_window_bounded(monkeypatch):
    gate = threading.Event()
    observed = {}

    class DB:
        def __init__(self):
            self.yielded = 0

        def features_of_type(self, _kind):
            for index in range(1000):
                self.yielded += 1
                yield SimpleNamespace(id=f"g{index}", index=index)

    class Factory:
        def __init__(self, _ctx):
            self.viable = True

        def open_root_scan(self):
            return db

        def close_root_scan(self, _db):
            return None

        def close(self):
            return None

    db = DB()

    def materialise(index, locus, _ctx):
        return MaterialisedLocus(index, locus, locus.id)

    def process(payload, _worker_ctx):
        if payload.submission_index == 0:
            assert gate.wait(3)
        return LocusResult(
            payload.submission_index, payload.locus_id,
            _Gene(payload.submission_index),
        )

    monkeypatch.setattr(locus_pipeline, "_ThreadLocalCtxFactory", Factory)
    monkeypatch.setattr(locus_pipeline, "materialise_locus", materialise)
    monkeypatch.setattr(
        locus_pipeline, "materialise_locus_with_factory", materialise,
    )
    monkeypatch.setattr(locus_pipeline, "process_locus_native", process)

    def release():
        time.sleep(0.1)
        observed["roots"] = db.yielded
        gate.set()

    releaser = threading.Thread(target=release)
    releaser.start()
    ctx = _ctx(db)
    count = parallel.parallel_step7(
        ["gene"], db, ctx, io.StringIO(),
        {"coding": {}, "non-coding": {}, "other": {}},
        threads=3, progress_every=0,
    )
    releaser.join()

    assert count == 1000
    assert observed["roots"] == 4
    assert ctx.args._step7_max_inflight_observed == 4


def test_nonreopenable_database_parent_materialises_then_fans_out(
        monkeypatch):
    loci = [SimpleNamespace(id=f"g{i}") for i in range(20)]
    main_thread = threading.current_thread().ident
    first_locus_gate = threading.Event()
    lock = threading.Lock()
    observed = {}
    worker_threads = set()

    class Connection:
        def close(self):
            return None

    class DB:
        dbfn = ":memory:"
        conn = Connection()

        def features_of_type(self, _kind):
            raise AssertionError("unsafe root scan used on shared connection")

        def features_of_type_stream_safe(self, _kind):
            for locus in loci:
                with lock:
                    observed["roots"] = observed.get("roots", 0) + 1
                yield locus

    db = DB()
    ctx = _ctx(db)
    ctx.args.step7_max_inflight = 4

    def materialise(index, locus, materialise_ctx):
        assert threading.current_thread().ident == main_thread
        assert materialise_ctx.l_feature_db is db
        with lock:
            observed["materialised"] = observed.get("materialised", 0) + 1
        return MaterialisedLocus(index, locus, locus.id)

    def process(payload, worker_ctx):
        with lock:
            worker_threads.add(threading.current_thread().ident)
        if payload.submission_index == 0:
            assert first_locus_gate.wait(3)
        delta = worker_ctx.state_journal.finish()
        return LocusResult(
            payload.submission_index, payload.locus_id,
            _Gene(payload.submission_index), delta=delta,
        )

    monkeypatch.setattr(locus_pipeline, "materialise_locus", materialise)
    monkeypatch.setattr(locus_pipeline, "process_locus_native", process)

    def release_first():
        time.sleep(0.1)
        with lock:
            observed["window_roots"] = observed.get("roots", 0)
            observed["window_materialised"] = observed.get("materialised", 0)
        first_locus_gate.set()

    releaser = threading.Thread(target=release_first)
    releaser.start()
    output = io.StringIO()

    count = parallel.parallel_step7(
        ["gene"], db, ctx, output,
        {"coding": {}, "non-coding": {}, "other": {}},
        threads=3, progress_every=0,
    )
    releaser.join()

    assert count == len(loci)
    assert output.getvalue() == "".join(f"{index}\n" for index in range(20))
    assert observed["window_roots"] == 4
    assert observed["window_materialised"] == 4
    assert ctx.args._step7_max_inflight_observed == 4
    assert ctx.args._step7_parallel_strategy == (
        "parallel_parent_materialise_worker_process"
    )
    assert "non-reopenable" in ctx.args._step7_materialisation_constraint
    assert not hasattr(ctx.args, "_step7_parallel_fallback_reason")
    assert main_thread not in worker_threads
    assert len(worker_threads) >= 2


def test_nonreopenable_parent_materialisation_failure_is_isolated(
        monkeypatch):
    loci = [SimpleNamespace(id=f"g{i}") for i in range(5)]

    class Connection:
        def close(self):
            return None

    class DB:
        dbfn = ":memory:"
        conn = Connection()

        def features_of_type_stream_safe(self, _kind):
            yield from loci

    db = DB()
    ctx = _ctx(db, window=3)

    def materialise(index, locus, _ctx):
        if index == 2:
            raise RuntimeError("injected materialisation failure")
        return MaterialisedLocus(index, locus, locus.id)

    def process(payload, worker_ctx):
        delta = worker_ctx.state_journal.finish()
        return LocusResult(
            payload.submission_index, payload.locus_id,
            _Gene(payload.submission_index), delta=delta,
        )

    monkeypatch.setattr(locus_pipeline, "materialise_locus", materialise)
    monkeypatch.setattr(locus_pipeline, "process_locus_native", process)
    output = io.StringIO()

    count = parallel.parallel_step7(
        ["gene"], db, ctx, output,
        {"coding": {}, "non-coding": {}, "other": {}},
        threads=3, progress_every=0,
    )

    assert count == len(loci)
    assert output.getvalue() == "0\n1\n3\n4\n"
    assert [record["locus_id"] for record in ctx.failure_records] == ["g2"]
    assert ctx.args._step7_max_inflight_observed <= 3


def test_root_scan_open_failure_closes_factory(monkeypatch):
    closed = []

    class DB:
        def features_of_type(self, _kind):
            return iter(())

    class Factory:
        viable = True

        def __init__(self, _ctx):
            pass

        def open_root_scan(self):
            raise RuntimeError("injected root scan failure")

        def close(self):
            closed.append(True)

    db = DB()
    monkeypatch.setattr(locus_pipeline, "_ThreadLocalCtxFactory", Factory)

    with pytest.raises(RuntimeError, match="root scan failure"):
        parallel.parallel_step7(
            ["gene"], db, _ctx(db), io.StringIO(),
            {"coding": {}, "non-coding": {}, "other": {}},
            threads=2, progress_every=0,
        )

    assert closed == [True]


def test_thread_local_factory_closes_all_reopened_handles(tmp_path):
    db_path = tmp_path / "features.db"
    db_path.write_bytes(b"placeholder")

    class Connection:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    class DB:
        dbfn = str(db_path)
        opened = []

        def __init__(self, _path=None):
            self.conn = Connection()
            if _path is not None:
                self.opened.append(self)

    parent = DB()
    ctx = _ctx(parent)
    ctx.ref_db = parent
    factory = locus_pipeline._ThreadLocalCtxFactory(ctx)
    scan = factory.open_root_scan()
    worker_ctx = factory.get()
    worker_handles = [worker_ctx.ref_db, worker_ctx.l_feature_db]
    factory.close()
    factory.close_root_scan(scan)

    assert scan is not parent
    assert all(handle.conn.closed for handle in worker_handles + [scan])
    assert parent.conn.closed is False


def test_thread_local_factory_closes_partial_reopens_on_failure(tmp_path):
    ref_path = tmp_path / "ref.db"
    liftoff_path = tmp_path / "liftoff.db"
    ref_path.touch()
    liftoff_path.touch()

    class Connection:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    class DB:
        opened = []

        def __init__(self, path):
            if os.fspath(path) == os.fspath(liftoff_path):
                raise RuntimeError("injected reopen failure")
            self.dbfn = os.fspath(path)
            self.conn = Connection()
            self.opened.append(self)

    def parent(path):
        db = object.__new__(DB)
        db.dbfn = os.fspath(path)
        db.conn = Connection()
        return db

    ctx = _ctx(parent(liftoff_path))
    ctx.ref_db = parent(ref_path)
    factory = locus_pipeline._ThreadLocalCtxFactory(ctx)

    with pytest.raises(RuntimeError, match="could not reopen"):
        factory.get()

    assert len(DB.opened) == 1
    assert DB.opened[0].conn.closed is True


def test_thread_local_factory_clones_inmemory_gffbase_connections():
    text = (
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\tx\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
    )
    parent = create_db(
        text, ":memory:", from_string=True, force=True, build_rtree=False,
    )
    ctx = _ctx(parent)
    ctx.ref_db = parent
    factory = locus_pipeline._ThreadLocalCtxFactory(ctx)

    assert factory.viable is True
    root_scan = factory.open_root_scan()
    worker_ctx = factory.get()
    assert root_scan is not parent
    assert worker_ctx.l_feature_db is not parent
    assert worker_ctx.ref_db is not parent
    assert [feature.id for feature in root_scan.features_of_type("gene")] == [
        "g1",
    ]
    assert [
        feature.id for feature in worker_ctx.l_feature_db.children(
            "g1", level=1, order_by="start",
        )
    ] == ["t1", "t2"]

    factory.close()
    factory.close_root_scan(root_scan)
    # Closing every duplicate must leave the owner usable.
    assert parent.count_features_of_type("gene") == 1


def test_thread_local_factory_clones_inmemory_gffbase_concurrently():
    records = ["##gff-version 3\n"]
    for index in range(32):
        start = index * 100 + 1
        end = start + 89
        records.extend([
            f"chr1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=g{index}\n",
            f"chr1\tx\tmRNA\t{start}\t{end}\t.\t+\t.\t"
            f"ID=t{index}b;Parent=g{index}\n",
            f"chr1\tx\tmRNA\t{start}\t{end}\t.\t+\t.\t"
            f"ID=t{index}a;Parent=g{index}\n",
        ])
    parent = create_db(
        "".join(records), ":memory:", from_string=True, force=True,
        build_rtree=False,
    )
    factory = locus_pipeline._ThreadLocalCtxFactory(_ctx(parent))
    worker_count = 8
    before_open = threading.Barrier(worker_count)
    before_query = threading.Barrier(worker_count)

    def inspect_gene_group(group):
        # Force connection duplication and queries to overlap across every
        # worker, rather than letting a tiny fixture reuse one pool thread.
        before_open.wait(timeout=10)
        worker_ctx = factory.get()
        before_query.wait(timeout=10)
        observed = []
        for gene_id in group:
            observed.append((
                gene_id,
                tuple(feature.id for feature in
                      worker_ctx.l_feature_db.children(
                          gene_id, level=1, order_by="start",
                      )),
            ))
        return (
            threading.get_ident(), id(worker_ctx.l_feature_db), observed,
        )

    groups = [
        [f"g{index}" for index in range(offset, 32, worker_count)]
        for offset in range(worker_count)
    ]
    try:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = [executor.submit(inspect_gene_group, group)
                       for group in groups]
            results = [future.result(timeout=20) for future in futures]
    finally:
        factory.close()

    assert len({thread_id for thread_id, _db_id, _rows in results}) == (
        worker_count
    )
    assert len({db_id for _thread_id, db_id, _rows in results}) == worker_count
    for _thread_id, _db_id, rows in results:
        for gene_id, child_ids in rows:
            index = gene_id.removeprefix("g")
            assert child_ids == (f"t{index}a", f"t{index}b")
    # Closing every worker duplicate must leave the owner catalog usable.
    assert parent.count_features_of_type("gene") == 32


def test_full_feature_batches_match_scalar_hierarchy():
    text = (
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1;Name=one\n"
        "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\texon\t1\t20\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tx\texon\t40\t80\t.\t+\t.\tID=e2;Parent=t1\n"
        "chr1\tx\tgene\t200\t300\t.\t-\t.\tID=g2\n"
    )
    db = create_db(
        text, ":memory:", from_string=True, force=True, build_rtree=False,
    )
    fetched = db.features_batched(["g2", "missing", "g1"])
    assert list(fetched) == ["g2", "g1"]
    assert fetched["g1"].attributes["Name"] == ["one"]

    batched = db.children_batched_features(
        ["g2", "g1"], order_by="start",
    )
    assert list(batched) == ["g2", "g1"]
    for anchor in ("g2", "g1"):
        scalar = tuple(db.children(anchor, order_by="start"))
        assert [str(feature) for feature in batched[anchor]] == [
            str(feature) for feature in scalar
        ]
        assert [feature.file_order for feature in batched[anchor]] == [
            feature.file_order for feature in scalar
        ]


def test_hierarchy_order_by_start_uses_gffutils_id_tiebreaker():
    text = (
        "##gff-version 3\n"
        "chr1\tmp\tmRNA\t1\t100\t.\t-\t.\tID=MP1\n"
        "chr1\tmp\tstop_codon\t10\t12\t.\t-\t0\t"
        "ID=stop;Parent=MP1\n"
        "chr1\tmp\tCDS\t10\t30\t.\t-\t2\tID=cds;Parent=MP1\n"
        "chr1\tmp\tCDS\t20\t40\t.\t-\t0\tID=later;Parent=MP1\n"
    )
    db = create_db(
        text, ":memory:", from_string=True, force=True, build_rtree=False,
    )

    ascending = list(db.children(
        "MP1", featuretype=("CDS", "stop_codon"), order_by="start",
    ))
    descending = list(db.children(
        "MP1", featuretype=("CDS", "stop_codon"), order_by="start",
        reverse=True,
    ))
    batched = db.children_batched_features(
        ["MP1"], featuretype=("CDS", "stop_codon"), order_by="start",
    )["MP1"]

    assert [feature.id for feature in ascending] == ["cds", "stop", "later"]
    assert [feature.id for feature in descending] == ["later", "cds", "stop"]
    assert [feature.id for feature in batched] == ["cds", "stop", "later"]


def test_gffbase_root_scan_survives_interleaved_child_queries(monkeypatch):
    monkeypatch.setattr(gffbase_interface, "_SCAN_BATCH_SIZE", 2)
    text = "##gff-version 3\n" + "".join(
        f"chr1\tx\tgene\t{i * 100 + 1}\t{i * 100 + 90}\t.\t+\t.\tID=g{i}\n"
        f"chr1\tx\tmRNA\t{i * 100 + 1}\t{i * 100 + 90}\t.\t+\t.\t"
        f"ID=t{i};Parent=g{i}\n"
        for i in range(3)
    )
    db = create_db(
        text, ":memory:", from_string=True, force=True, build_rtree=False,
    )

    observed = []
    for gene in db.features_of_type("gene"):
        observed.append(gene.id)
        assert len(list(db.children(gene, level=1))) == 1

    assert observed == ["g0", "g1", "g2"]


def test_chunk_decoder_handles_split_utf8_records_without_blob():
    decoder = GFF3ChunkDecoder(max_chunk_bytes=7)
    records = []
    for offset in range(0, len(GFF), 5):
        records.extend(decoder.feed(GFF[offset:offset + 5]))
    records.extend(decoder.finish())

    assert [record.featuretype for record in records] == ["mRNA", "CDS"]
    assert decoder.byte_count == len(GFF)
    assert decoder.feature_count == 2
    assert decoder.pending_bytes == 0
    assert decoder.directives == ["##gff-version 3"]


def test_record_stream_closes_connection_on_ingest_failure(monkeypatch):
    class Connection:
        closed = False

        def close(self):
            self.closed = True

    connection = Connection()

    def fail(*_args, _connection_out, **_kwargs):
        _connection_out.append(connection)
        raise RuntimeError("injected ingest failure")

    monkeypatch.setattr(gffbase_ingest, "_from_record_stream_impl", fail)

    with pytest.raises(RuntimeError, match="ingest failure"):
        gffbase_ingest.from_record_stream(())

    assert connection.closed is True


class _Process:
    def __init__(self, stdout=GFF, stderr=b"", returncode=0):
        self.stdout = io.BytesIO(stdout)
        self.stderr = io.BytesIO(stderr)
        self._returncode = returncode
        self.returncode = None

    def wait(self, timeout=None):
        self.returncode = self._returncode
        return self._returncode

    def poll(self):
        return self.returncode

    def terminate(self):
        self.returncode = -15

    def kill(self):
        self.returncode = -9


def test_miniprot_stream_publishes_queryable_database(tmp_path, monkeypatch):
    fsynced = []
    monkeypatch.setattr(
        run_miniprot.subprocess, "Popen", lambda *args, **kwargs: _Process(),
    )
    monkeypatch.setattr(
        "lifton.output_transaction._fsync_directory",
        lambda path: fsynced.append(path),
    )
    artifact = run_miniprot.run_miniprot_streaming_db(
        os.fspath(tmp_path) + "/",
        SimpleNamespace(mp_options="", threads=2), "target.fa", "ref.fa",
        chunk_size=11,
    )

    assert artifact.byte_count == len(GFF)
    assert artifact.feature_count == 2
    assert artifact.returncode == 0
    assert fsynced == [Path(artifact.database_path).parent]
    assert os.path.isfile(artifact.database_path)
    assert not os.path.exists(artifact.database_path + ".partial")
    from lifton.gffbase_adapter import open_database_path
    opened = open_database_path(artifact.database_path)
    assert opened.count_features_of_type("mRNA") == 1
    opened.conn.close()


def test_miniprot_stream_failure_preserves_old_database(tmp_path, monkeypatch):
    destination = tmp_path / "miniprot" / "miniprot.duckdb"
    destination.parent.mkdir()
    destination.write_bytes(b"old-result")
    partial_wal = Path(os.fspath(destination) + ".partial.wal")
    partial_wal.write_bytes(b"stale-wal")
    monkeypatch.setattr(
        run_miniprot.subprocess, "Popen",
        lambda *args, **kwargs: _Process(stderr=b"prefix ERROR suffix"),
    )

    with pytest.raises(RuntimeError, match="ERROR"):
        run_miniprot.run_miniprot_streaming_db(
            os.fspath(tmp_path) + "/",
            SimpleNamespace(mp_options="", threads=1),
            "target.fa", "ref.fa", chunk_size=9,
        )

    assert destination.read_bytes() == b"old-result"
    assert not os.path.exists(os.fspath(destination) + ".partial")
    assert not partial_wal.exists()


def test_miniprot_stream_refuses_stale_destination_wal(tmp_path, monkeypatch):
    destination = tmp_path / "miniprot" / "miniprot.duckdb"
    destination.parent.mkdir()
    destination.write_bytes(b"old-result")
    destination_wal = Path(os.fspath(destination) + ".wal")
    destination_wal.write_bytes(b"old-wal")

    def unexpected_popen(*_args, **_kwargs):
        raise AssertionError("miniprot must not start while a stale WAL exists")

    monkeypatch.setattr(run_miniprot.subprocess, "Popen", unexpected_popen)

    with pytest.raises(RuntimeError, match="destination WAL exists"):
        run_miniprot.run_miniprot_streaming_db(
            os.fspath(tmp_path) + "/",
            SimpleNamespace(mp_options="", threads=1),
            "target.fa", "ref.fa", chunk_size=9,
        )

    assert destination.read_bytes() == b"old-result"
    assert destination_wal.read_bytes() == b"old-wal"
    assert not Path(os.fspath(destination) + ".partial").exists()


def test_miniprot_facade_align_all_into_streams_chunks(monkeypatch):
    process = _Process()
    monkeypatch.setattr(
        miniprot_facade.subprocess, "Popen",
        lambda *args, **kwargs: process,
    )
    index = MiniprotIndex("target.fa", ref_proteins_path="ref.fa")
    chunks = []
    count = index.align_all_into(chunks.append, chunk_size=13)

    assert count == len(GFF)
    assert b"".join(chunks) == GFF
    assert max(map(len, chunks)) <= 13
    assert index._cached_raw_bytes is None
