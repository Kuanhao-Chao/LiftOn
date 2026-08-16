from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.compare import release_rescore


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _fixture(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    input_path = source / "reference.fa"
    input_path.write_text(">chr1\nAAAA\n", encoding="utf-8")
    input_record = release_rescore._artifact_record(input_path)
    versions = {}
    output_records = {}
    for label in ("candidate", "reference"):
        output = source / f"{label}.gff3"
        output.write_text(
            f"##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID={label}\n",
            encoding="utf-8",
        )
        fingerprints = {
            "byte_sha256": release_rescore.release_evaluation.sha256_file(output),
            "semantic_sha256": label[0] * 64,
        }
        versions[label] = {
            "output_gff": str(output.resolve()),
            "fingerprints": fingerprints,
            "summary": {
                "transcripts_tsv": f"sealed/{label}.tsv",
                "target_truth": {"source": "sealed"},
            },
            "evaluation_artifacts": {},
        }
        output_records[label] = {
            "path": str(output.resolve()),
            "size": output.stat().st_size,
            "fingerprints": fingerprints,
        }
    pair = {
        "panel": "subset",
        "benchmark": "demo",
        "repetition": 1,
        "inputs": {"ref_fa": input_record},
        "versions": versions,
    }
    pair_path = source / "pair_result.json"
    _write_json(pair_path, pair)
    pair_record = release_rescore._artifact_record(pair_path)

    root = tmp_path / "overlay"
    record_dir = root / "records" / "subset" / "demo" / "repetition-01"
    overlay_versions = {}
    for label in ("candidate", "reference"):
        transcript = record_dir / "evaluation" / f"{label}.transcripts.tsv"
        transcript.parent.mkdir(parents=True, exist_ok=True)
        transcript.write_text(f"label\n{label}\n", encoding="utf-8")
        overlay_versions[label] = {
            "summary": {
                "transcripts_tsv": f"evaluation/{label}.transcripts.tsv",
                "rescored": label,
            },
            "transcripts_tsv": {
                "path": f"evaluation/{label}.transcripts.tsv",
                "size": transcript.stat().st_size,
                "sha256": release_rescore.release_evaluation.sha256_file(transcript),
            },
        }
    tooling_file = tmp_path / "evaluator.py"
    tooling_file.write_text("# evaluator snapshot\n", encoding="utf-8")
    tooling = {"evaluator": release_rescore._artifact_record(tooling_file)}
    tooling_sha256 = release_rescore._sha256_text(
        release_rescore._canonical_json(tooling)
    )
    record = {
        "schema_version": release_rescore.OVERLAY_SCHEMA_VERSION,
        "kind": release_rescore.RECORD_KIND,
        "key": {"panel": "subset", "benchmark": "demo", "repetition": 1},
        "source": {
            "pair_result": pair_record,
            "inputs": pair["inputs"],
            "outputs": output_records,
        },
        "tooling_sha256": tooling_sha256,
        "versions": overlay_versions,
    }
    record_path = record_dir / "overlay.json"
    _write_json(record_path, record)
    manifest = {
        "schema_version": release_rescore.OVERLAY_SCHEMA_VERSION,
        "kind": release_rescore.OVERLAY_KIND,
        "created_at": "2026-08-10T00:00:00Z",
        "run_roots": [str(source.resolve())],
        "source_pair_count": 1,
        "tooling": tooling,
        "tooling_sha256": tooling_sha256,
        "runtime": {"python": "test"},
        "git_state": {"head": "test"},
        "records": [{
            "key": dict(record["key"]),
            "source_pair": pair_record,
            "record": str(record_path.relative_to(root)),
            "record_sha256": release_rescore.release_evaluation.sha256_file(
                record_path
            ),
        }],
    }
    _write_json(root / "manifest.json", manifest)
    return pair_path, pair, root, record_dir


def test_apply_overlay_replaces_only_evaluation_evidence(tmp_path):
    pair_path, pair, root, _record_dir = _fixture(tmp_path)

    transformed, evidence = release_rescore.apply_overlays(
        [(pair_path, pair)], [root],
    )

    rescored = transformed[0][1]
    for label in ("candidate", "reference"):
        version = rescored["versions"][label]
        assert version["summary"]["rescored"] == label
        assert version["summary"]["target_truth"] == {"source": "sealed"}
        assert Path(version["summary"]["transcripts_tsv"]).is_absolute()
        assert version["output_gff"] == pair["versions"][label]["output_gff"]
    assert evidence["validated"] is True
    assert evidence["pair_count"] == 1


def test_overlay_rejects_modified_transcript_artifact(tmp_path):
    pair_path, pair, root, record_dir = _fixture(tmp_path)
    transcript = record_dir / "evaluation" / "candidate.transcripts.tsv"
    transcript.write_text("tampered\n", encoding="utf-8")

    with pytest.raises(ValueError, match="artifact was modified"):
        release_rescore.apply_overlays([(pair_path, pair)], [root])


def test_overlay_requires_exact_pair_coverage(tmp_path):
    _pair_path, _pair, root, _record_dir = _fixture(tmp_path)

    with pytest.raises(ValueError, match="coverage does not match"):
        release_rescore.apply_overlays([], [root])


def test_write_overlay_rescores_mrna_and_transcript_outputs(tmp_path):
    run_root = tmp_path / "runs"
    pair_dir = run_root / "demo"
    pair_dir.mkdir(parents=True)
    sequence = "ATG" + "AAA" * 8 + "TAA"
    ref_fa = tmp_path / "reference.fa"
    tgt_fa = tmp_path / "target.fa"
    for path in (ref_fa, tgt_fa):
        path.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
    ref_gff = tmp_path / "reference.gff3"
    ref_gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\ttranscript\t1\t30\t.\t+\t.\tID=tx1;Parent=g1\n"
        "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=ex1;Parent=tx1\n"
        "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1\n",
        encoding="utf-8",
    )
    versions = {}
    for label, transcript_type in (
        ("candidate", "mRNA"),
        ("reference", "transcript"),
    ):
        output = tmp_path / f"{label}.gff3"
        output.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
            f"chr1\ttest\t{transcript_type}\t1\t30\t.\t+\t.\t"
            "ID=tx1;Parent=g1\n"
            "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=ex1;Parent=tx1\n"
            "chr1\ttest\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1\n",
            encoding="utf-8",
        )
        versions[label] = {
            "output_gff": str(output.resolve()),
            "fingerprints": release_rescore.release_evaluation.gff3_fingerprints(
                output
            ),
            "profile": {},
            "summary": {"species": "synthetic", "cross_species": False},
        }
    pair = {
        "schema_version": release_rescore.release_evaluation.SCHEMA_VERSION,
        "panel": "subset",
        "benchmark": "demo",
        "repetition": 1,
        "inputs": {
            name: release_rescore._artifact_record(path)
            for name, path in {
                "ref_fa": ref_fa,
                "ref_gff": ref_gff,
                "tgt_fa": tgt_fa,
            }.items()
        },
        "versions": versions,
    }
    pair_path = pair_dir / "pair_result.json"
    _write_json(pair_path, pair)
    overlay_root = tmp_path / "overlay-output"

    release_rescore.write_overlay(
        [run_root], overlay_root, threads=1, allow_incomplete=True,
    )
    transformed, _evidence = release_rescore.apply_overlays(
        [(pair_path, pair)], [overlay_root],
    )

    for label in ("candidate", "reference"):
        summary = transformed[0][1]["versions"][label]["summary"]
        assert summary["n_reference_coding"] == 1
        assert summary["n_recovered_coding"] == 1
        assert summary["completeness_coding"] == 1.0


def test_tooling_records_cover_the_first_party_scoring_tree():
    records = release_rescore._tooling_records()
    package_root = Path(release_rescore.build_controller.REPO_ROOT) / "lifton"
    expected = {
        f"lifton/{path.relative_to(package_root).as_posix()}"
        for path in package_root.rglob("*.py")
        if path.stat().st_size > 0
    }

    assert expected <= set(records)
    assert {
        "biopython", "gffutils", "intervaltree", "numpy", "parasail",
        "pyfaidx", "pysam",
    } == set(release_rescore._runtime_record()["packages"])
