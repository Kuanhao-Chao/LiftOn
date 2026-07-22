from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from benchmarks.compare import (
    build_controller,
    campaign_profiles,
    protocol_analysis,
    release_report,
)


def test_legacy_profile_is_exactly_the_frozen_34_17_5_campaign():
    profile = campaign_profiles.load_profile("canonical-v1")
    counts = campaign_profiles.campaign_counts(profile)
    specification = campaign_profiles.campaign_spec(
        profile,
        legacy_v1=True,
    )

    assert specification == release_report.canonical_campaign_spec()
    assert {
        panel: len(row["ids"])
        for panel, row in specification["panels"].items()
    } == {"subset": 34, "full": 17, "e2e": 5}
    assert counts["cells"] == 224
    assert counts["arms"] == 448
    assert "t1_soybean_w82_to_lee" not in specification["panels"]["subset"]["ids"]


def test_v2_profile_has_explicit_release_truth_and_shared_protocol_matrices():
    profile = campaign_profiles.load_profile("canonical-v2")
    counts = campaign_profiles.campaign_counts(profile)
    specification = campaign_profiles.campaign_spec(profile)
    by_id = {row["id"]: row for row in counts["cases"]}

    assert {
        panel: len(row["ids"])
        for panel, row in specification["panels"].items()
    } == {"subset": 40, "full": 22, "e2e": 11}
    assert len(specification["matrix"]) == 7
    assert counts["campaigns"] == 16
    assert counts["cells"] == 411
    assert counts["arms"] == 751
    assert by_id["synthetic-truth"]["cells"] == 8
    assert by_id["subset-deep-diagnostic"]["cells"] == 8
    assert by_id["full-deep-diagnostic"]["cells"] == 8
    assert {
        item
        for row in specification["matrix"]
        for item in row["ids"]
        if item.startswith("v2_deep_")
    } == set()
    assert by_id["v2_protocol_thread_scaling_bee"]["cells"] == len(
        protocol_analysis.protocol_cases("thread_scaling")
    ) == 36
    assert by_id["v2_protocol_io_modes_arabidopsis"]["cells"] == len(
        protocol_analysis.protocol_cases("io_modes")
    ) == 32
    assert release_report._normalize_campaign_spec(
        specification
    ) == specification


def test_profile_loader_rejects_unknown_fields_and_unbalanced_release_runs(
        tmp_path):
    source = json.loads(
        campaign_profiles.DEFAULT_PROFILE_REGISTRY.read_text()
    )
    source["profiles"][0]["campaigns"][0]["surprise"] = True
    path = tmp_path / "unknown.json"
    path.write_text(json.dumps(source))

    with pytest.raises(ValueError, match="unknown=.*surprise"):
        campaign_profiles.load_registry(path)

    del source["profiles"][0]["campaigns"][0]["surprise"]
    source["profiles"][0]["campaigns"][0]["repetitions"] = 3
    path.write_text(json.dumps(source))
    with pytest.raises(ValueError, match="even AB/BA"):
        campaign_profiles.load_registry(path)


def test_profile_campaign_spec_detects_matrix_panel_disagreement():
    specification = release_report.canonical_campaign_spec(
        profile_id="canonical-v2",
    )
    specification["panels"]["subset"]["ids"].pop()

    with pytest.raises(ValueError, match="do not match"):
        release_report._normalize_campaign_spec(specification)


def test_controller_profile_flags_fill_case_contract_and_reject_conflicts(
        tmp_path):
    parser = build_controller.build_parser()
    args = parser.parse_args([
        "start",
        "--campaign-profile", "canonical-v1",
        "--campaign-id", "subset",
        "--candidate-root", str(tmp_path / "candidate"),
        "--candidate-sha", "a" * 40,
        "--reference-root", str(tmp_path / "reference"),
        "--reference-sha", "b" * 40,
        "--paired-lifton-executable", sys.executable,
        "--dry-run",
    ])
    identity = build_controller._campaign_case_from_args(args)

    assert identity["profile_id"] == "canonical-v1"
    assert identity["case"]["id"] == "subset"
    assert args.stage == "paired-subset"
    assert len(args.ids) == 34
    assert args.paired_repetitions == 4
    assert args.threads_per_cell == 8

    conflicting = parser.parse_args([
        "start",
        "--campaign-profile", "canonical-v1",
        "--campaign-id", "subset",
        "--stage", "paired-full",
        "--dry-run",
    ])
    with pytest.raises(ValueError, match="conflict"):
        build_controller._campaign_case_from_args(conflicting)


def test_profile_registry_requires_profile_and_case_together():
    parser = build_controller.build_parser()
    args = parser.parse_args([
        "start",
        "--campaign-profile", "canonical-v1",
    ])
    with pytest.raises(ValueError, match="both"):
        build_controller._campaign_case_from_args(args)
