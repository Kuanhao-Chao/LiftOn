"""Tests for the generated canonical-v2 coverage and future gap matrix."""
from benchmarks.compare import canonical_v2_coverage


def test_coverage_matrix_matches_the_frozen_campaign_shape():
    document = canonical_v2_coverage.build_coverage()

    assert document["campaign_summary"]["campaigns"] == 16
    assert document["campaign_summary"]["cells"] == 411
    assert document["campaign_summary"]["arms"] == 751
    assert document["release_matrix"] == {
        "subset_ids": 40,
        "full_ids": 22,
        "e2e_ids": 11,
        "synthetic_ids": 2,
        "repetitions": [4],
    }
    assert len(document["scenarios"]) == 12
    assert document["future_backlog"]["items"]


def test_committed_coverage_outputs_are_current():
    document = canonical_v2_coverage.build_coverage()

    assert canonical_v2_coverage.DEFAULT_JSON.read_text() == (
        canonical_v2_coverage.render_json(document)
    )
    assert canonical_v2_coverage.DEFAULT_MARKDOWN.read_text() == (
        canonical_v2_coverage.render_markdown(document)
    )
