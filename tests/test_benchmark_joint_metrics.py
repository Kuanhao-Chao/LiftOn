import csv

from benchmarks.compare import joint_metrics


FIELDS = [
    "ref_mrna_id", "recovered", "is_coding", "protein_identity",
    "intron_chain_exact", "exon_sn", "exon_sp", "orf_valid",
]


def _write_rows(tmp_path, tool, rows):
    path = tmp_path / f"{tool}.transcripts.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _row(feature_id, identity, **overrides):
    row = {
        "ref_mrna_id": feature_id,
        "recovered": "1",
        "is_coding": "1",
        "protein_identity": str(identity),
        "intron_chain_exact": "1",
        "exon_sn": "0.8",
        "exon_sp": "0.9",
        "orf_valid": "1",
    }
    row.update({key: str(value) for key, value in overrides.items()})
    return row


def test_joint_metrics_include_thresholds_bootstrap_and_structure(tmp_path):
    _write_rows(tmp_path, "lifton_devel", [
        _row("a", 0.96),
        _row("b", 0.80, intron_chain_exact=0, exon_sn=0.6, orf_valid=0),
    ])
    _write_rows(tmp_path, "liftoff", [
        _row("a", 0.90),
        _row("b", 0.70),
    ])

    result = joint_metrics.compute_joint_metrics(str(tmp_path), 4)

    assert result["recall_at_0.5"]["lifton_devel"] == 0.5
    assert result["recall_at_0.75"]["lifton_devel"] == 0.5
    assert result["recall_at_0.9"]["lifton_devel"] == 0.25
    assert result["recall_at_0.95"]["lifton_devel"] == 0.25
    paired = result["devel_vs_liftoff_common"]
    assert paired["meanpi_delta"] == 0.08
    assert paired["n_common"] == 2
    assert paired["bootstrap_95ci"] == joint_metrics._paired(
        {"a": 0.96, "b": 0.80}, {"a": 0.90, "b": 0.70},
    )["bootstrap_95ci"]
    assert result["structural"]["lifton_devel"]["intron_chain_exact_recall"] == 0.25
    assert result["structural"]["lifton_devel"]["orf_valid_recall"] == 0.25
    assert result["structural"]["lifton_devel"]["exon_sn_mean"] == 0.7


def test_legacy_tsv_without_structural_columns_still_scores(tmp_path):
    path = tmp_path / "lifton_devel.transcripts.tsv"
    path.write_text(
        "ref_mrna_id\trecovered\tis_coding\tprotein_identity\n"
        "a\t1\t1\t0.9\n"
    )

    result = joint_metrics.compute_joint_metrics(str(tmp_path), 2)

    assert result["covpi"]["lifton_devel"] == 0.45
    assert "structural" not in result
