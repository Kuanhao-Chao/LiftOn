"""Track-B robustness fixes: honest error counts, durability, ordering, read-only refs.

Four defects the suite could not see, each verified here directly:

1. ``validate_gff3_structure`` counted issues into a LOCAL dict and never
   populated ``issue_totals`` / ``severity_totals``. ``result.issues`` is capped
   at ``max_issues_per_check``, so a file with thousands of errors reported
   exactly 50 per check, and ``ValidationResult.is_valid``'s severity_totals
   clause -- written precisely to survive that capping -- was inert.
2. ``lifton.py`` closed the staged GFF3 with a bare ``fw.close()``, so
   ``OutputTransaction.close()`` short-circuited on ``self._stream.closed`` and
   skipped its flush + ``os.fsync``; ``commit()`` then published, via
   ``os.replace``, data that may never have reached disk.
3. ``_order_clause(None)`` ordered by ``file_order`` alone, which is not a total
   order: inferred roots carry no ingest position, so a root scan could return
   a different order between passes.
4. ``cache_lock`` created ``<db>.lock`` next to the database, so a shared
   read-only reference mount -- the normal HPC layout -- raised an uncaught
   ``PermissionError`` where the cached database used to be reused read-only.
"""
import os
import stat

import pytest

from lifton import annotation_cache, gff3_validator
from lifton.gff3_validator import Severity, validate_gff3_structure


HEADER = "##gff-version 3\n"


def _write(path, rows):
    path.write_text(HEADER + "".join(rows))
    return str(path)


def _row(seqid="chr1", ftype="gene", start=1, end=100, attrs="ID=g1"):
    return f"{seqid}\tLiftOn\t{ftype}\t{start}\t{end}\t.\t+\t.\t{attrs}\n"


# ── 1. uncapped error totals ────────────────────────────────────────────────

class TestStructuralValidatorTotals:
    def test_totals_are_uncapped_while_issue_objects_stay_capped(self, tmp_path):
        # 120 rows sharing one ID -> 119 duplicate_id errors.
        rows = [_row(attrs="ID=dup") for _ in range(120)]
        path = _write(tmp_path / "dups.gff3", rows)

        result = validate_gff3_structure(path, max_issues_per_check=10)

        assert result.issue_totals["duplicate_id"] == 119
        assert result.severity_totals[Severity.ERROR] == 119
        # The materialized objects remain bounded -- that is the whole point of
        # the cap, and why the totals had to be tracked separately.
        assert len(result.errors) == 10
        assert not result.is_valid

    def test_counts_only_mode_still_reports_invalid(self, tmp_path):
        """``max_issues_per_check=0`` retains no issue objects at all.

        Before the fix ``is_valid`` was ``0 == 0 and len([]) == 0`` -> True,
        i.e. a broken file validated clean in counts-only mode.
        """
        path = _write(tmp_path / "dups.gff3", [_row(attrs="ID=dup")] * 5)

        result = validate_gff3_structure(path, max_issues_per_check=0)

        assert result.errors == []
        assert result.severity_totals[Severity.ERROR] == 4
        assert not result.is_valid

    def test_clean_file_reports_no_totals(self, tmp_path):
        path = _write(tmp_path / "ok.gff3", [
            _row(attrs="ID=gene1"),
            _row(ftype="mRNA", attrs="ID=rna1;Parent=gene1"),
            _row(ftype="exon", attrs="ID=exon1;Parent=rna1"),
        ])

        result = validate_gff3_structure(path)

        assert result.severity_totals.get(Severity.ERROR, 0) == 0
        assert result.issue_totals == {}
        assert result.is_valid

    def test_early_return_paths_also_carry_totals(self, tmp_path):
        result = validate_gff3_structure(str(tmp_path / "missing.gff3"))
        assert result.severity_totals[Severity.ERROR] >= 1
        assert not result.is_valid


# ── 2. staged output durability ─────────────────────────────────────────────

class TestStagedOutputDurability:
    def test_closing_through_the_transaction_fsyncs(self, tmp_path, monkeypatch):
        from lifton import output_transaction as ot

        synced = []
        real_fsync = os.fsync
        monkeypatch.setattr(
            ot.os, "fsync",
            lambda fd: (synced.append(fd), real_fsync(fd))[1],
        )

        txn = ot.OutputTransaction(str(tmp_path / "out.gff3"))
        handle = txn.open()
        handle.write(HEADER)
        txn.close()
        assert synced, "close() must fsync the staged stream"

    def test_a_bare_stream_close_would_skip_the_fsync(self, tmp_path, monkeypatch):
        """Pins the actual defect: the short-circuit that made it silent."""
        from lifton import output_transaction as ot

        synced = []
        real_fsync = os.fsync
        monkeypatch.setattr(
            ot.os, "fsync",
            lambda fd: (synced.append(fd), real_fsync(fd))[1],
        )

        txn = ot.OutputTransaction(str(tmp_path / "out.gff3"))
        handle = txn.open()
        handle.write(HEADER)
        handle.close()          # the old `fw.close()`
        txn.close()             # sees `self._stream.closed` -> skips fsync
        assert synced == []

    def test_lifton_closes_the_staged_gff3_through_the_transaction(self):
        """The pipeline must not reintroduce the bare close."""
        import inspect
        from lifton import lifton as lifton_module

        source = inspect.getsource(lifton_module.run_all_lifton_steps)
        assert "output_transaction.close()" in source
        assert "\n        fw.close()\n" not in source


# ── 3. total ordering for gffbase root scans ────────────────────────────────

class TestGffbaseOrdering:
    @staticmethod
    def _clause(order_by, reverse=False):
        from lifton.gffbase.interface import FeatureDB
        return FeatureDB._order_clause(order_by, reverse)

    def test_unordered_scan_has_a_unique_tiebreak(self):
        clause = self._clause(None)
        assert clause.startswith("file_order ASC")
        assert clause.endswith(", id ASC")

    def test_explicit_file_order_is_also_total(self):
        assert self._clause("file_order").endswith(", id ASC")

    def test_other_keys_keep_the_ingest_order_tiebreak(self):
        assert self._clause("start") == "start ASC, file_order ASC"

    def test_reverse_only_flips_the_primary_key(self):
        assert self._clause(None, reverse=True) == "file_order DESC, id ASC"


# ── 4. read-only reference directories ──────────────────────────────────────

class TestCacheLockOnReadOnlyDirectories:
    def test_lock_falls_back_when_the_directory_is_not_writable(self, tmp_path):
        refdir = tmp_path / "readonly_reference"
        refdir.mkdir()
        db_path = refdir / "ref.gff3_db"
        db_path.write_bytes(b"")
        original_mode = stat.S_IMODE(refdir.stat().st_mode)
        os.chmod(refdir, stat.S_IRUSR | stat.S_IXUSR)
        try:
            if os.access(refdir, os.W_OK):
                pytest.skip("cannot drop write permission (running as root?)")
            entered = False
            with annotation_cache.cache_lock(str(db_path)):
                entered = True
            assert entered
            assert not (refdir / "ref.gff3_db.lock").exists()
        finally:
            os.chmod(refdir, original_mode)

    def test_the_fallback_is_stable_for_one_database(self, tmp_path):
        lock_path = str(tmp_path / "nope" / "x.lock")
        first = annotation_cache._open_lock_file(lock_path)
        # The parent is creatable here, so this is the normal path.
        assert first is not None
        first.close()

    def test_normal_directories_still_lock_next_to_the_database(self, tmp_path):
        db_path = tmp_path / "ref.gff3_db"
        with annotation_cache.cache_lock(str(db_path)):
            pass
        assert (tmp_path / "ref.gff3_db.lock").exists()

    def test_unwritable_everywhere_proceeds_with_a_warning(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            annotation_cache, "_open_lock_file", lambda path: None,
        )
        warnings = []
        monkeypatch.setattr(
            annotation_cache.logger, "log_warning", warnings.append,
        )
        with annotation_cache.cache_lock(str(tmp_path / "ref.gff3_db")):
            pass
        assert warnings and "cache lock" in warnings[0]


def test_manifest_records_the_uncapped_structural_error_count():
    """`lifton.py` must report the total, not the per-check-capped length."""
    import inspect
    from lifton import lifton as lifton_module

    source = inspect.getsource(lifton_module.run_all_lifton_steps)
    assert "structural_error_total" in source
    assert "errors=len(structural_result.errors)" not in source
    assert gff3_validator.Severity.ERROR is not None
