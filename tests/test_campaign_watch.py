import json
from pathlib import Path

from benchmarks.compare import campaign_watch


def test_snapshot_aggregates_active_cells_and_resource_alerts(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    cell_dir = run_dir / "cells" / "cell-a"
    cell_dir.mkdir(parents=True)
    (cell_dir / "status.json").write_text(json.dumps({
        "state": "running",
        "attempts": 1,
        "started_ns": 1,
        "updated_at": "2026-08-10T00:00:00Z",
    }))
    plan = {
        "run_id": "run",
        "stage": "paired-full",
        "fingerprint": "f" * 64,
        "policy": {"max_active": 6},
        "cells": [{
            "id": "cell-a",
            "benchmark": "bee",
            "repetition": 1,
            "threads": 8,
            "full_job": True,
            "cell_dir": str(cell_dir),
        }],
    }
    monkeypatch.setattr(campaign_watch.build_controller, "load_plan", lambda _: plan)
    monkeypatch.setattr(
        campaign_watch.build_controller,
        "cell_session_name",
        lambda *_: "cell-session",
    )
    monkeypatch.setattr(
        campaign_watch.build_controller,
        "tmux_has_session",
        lambda _: True,
    )
    monkeypatch.setattr(
        campaign_watch.build_controller,
        "_worker_identity_matches",
        lambda _: False,
    )
    monkeypatch.setattr(
        campaign_watch,
        "_resource_snapshot",
        lambda: {"load1": 120.0, "available_gib": 300.0},
    )
    monkeypatch.setattr(campaign_watch, "_free_space_gib", lambda _: 10.0)

    result = campaign_watch.snapshot(
        [run_dir],
        max_active=0,
        max_worker_threads=4,
        load1_limit=112.0,
        min_available_gib=384.0,
    )

    assert result["totals"] == {
        "active_cells": 1,
        "active_threads": 8,
        "counts": {"running": 1},
    }
    assert len(result["alerts"]) == 4
    assert result["runs"][0]["active"][0]["session_alive"] is True


def test_write_snapshot_replaces_latest_and_appends_history(tmp_path, monkeypatch):
    monkeypatch.setattr(
        campaign_watch,
        "snapshot",
        lambda *_args, **_kwargs: {"schema_version": 1, "alerts": []},
    )
    latest = tmp_path / "status" / "latest.json"
    history = tmp_path / "status" / "history.jsonl"

    campaign_watch.write_snapshot(
        [Path("run")], latest=latest, history=history,
    )

    assert json.loads(latest.read_text())["schema_version"] == 1
    assert len(history.read_text().splitlines()) == 1


def test_snapshot_uses_persisted_session_and_tmux_normalized_name(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    cell_dir = run_dir / "cells" / "cell-a"
    cell_dir.mkdir(parents=True)
    (cell_dir / "status.json").write_text(json.dumps({
        "state": "running",
        "session": "live.session",
        "worker_pid": 123,
        "worker_start_ticks": 456,
    }))
    plan = {
        "run_id": "run",
        "stage": "paired-full",
        "fingerprint": "f" * 64,
        "policy": {},
        "cells": [{
            "id": "cell-a",
            "benchmark": "bee",
            "repetition": 1,
            "threads": 8,
            "full_job": True,
            "cell_dir": str(cell_dir),
        }],
    }
    calls = []
    monkeypatch.setattr(campaign_watch.build_controller, "load_plan", lambda _: plan)
    monkeypatch.setattr(
        campaign_watch.build_controller,
        "tmux_has_session",
        lambda name: calls.append(name) or name == "live_session",
    )
    monkeypatch.setattr(
        campaign_watch.build_controller,
        "_worker_identity_matches",
        lambda _: True,
    )
    monkeypatch.setattr(
        campaign_watch,
        "_resource_snapshot",
        lambda: {"load1": 1.0, "available_gib": 500.0},
    )
    monkeypatch.setattr(campaign_watch, "_free_space_gib", lambda _: 10.0)

    result = campaign_watch.snapshot([run_dir])
    active = result["runs"][0]["active"][0]

    assert active["tmux_session"] == "live.session"
    assert active["session_alive"] is True
    assert calls == ["live.session", "live_session"]


def test_write_snapshot_records_evidence_manifest(tmp_path, monkeypatch):
    payload = {
        "schema_version": 1,
        "observed_at": "2026-08-10T00:00:00Z",
        "host": "host",
        "limits": {"max_active": 16},
        "totals": {"active_cells": 1},
        "alerts": [],
        "runs": [{
            "run_id": "run",
            "run_dir": str(tmp_path / "run"),
            "stage": "subset",
            "plan_fingerprint": "f" * 64,
            "policy": {"max_active": 16},
        }],
    }
    monkeypatch.setattr(campaign_watch, "snapshot", lambda *_args, **_kwargs: payload)
    latest = tmp_path / "latest.json"
    evidence = tmp_path / "evidence.json"

    campaign_watch.write_snapshot(
        [Path("run")], latest=latest, evidence=evidence,
    )

    document = json.loads(evidence.read_text())
    assert document["kind"] == "campaign_evidence_manifest"
    assert document["run_plans"][0]["plan_fingerprint"] == "f" * 64
    assert json.loads(latest.read_text())["evidence_manifest"] == str(
        evidence.resolve()
    )
def test_cpu_iowait_uses_the_fifth_proc_stat_counter():
    stat = "cpu  100 0 50 800 50 0 0 0 20 5\n"

    assert campaign_watch._cpu_iowait_percent_from_stat(stat) == 5.0
