from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.compare import release_evaluation as release


def _inputs(tmp_path: Path, panel: str) -> release.PanelInputs:
    tmp_path.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name in (
        "ref_fa",
        "ref_gff",
        "tgt_fa",
        "liftoff_gff",
        "miniprot_gff",
        "transcripts_fa",
        "proteins_fa",
    ):
        path = tmp_path / name
        path.write_text(f"{name}\n")
        paths[name] = path
    if panel == "subset":
        paths["transcripts_fa"] = None
        paths["proteins_fa"] = None
    return release.PanelInputs(
        benchmark="demo",
        panel=panel,
        species="Example",
        cross_species=False,
        annotation_database="RefSeq",
        **paths,
    )


def _source(tmp_path: Path) -> release.SourceSpec:
    return release.SourceSpec(
        label="candidate",
        root=tmp_path,
        sha="a" * 40,
        lifton_executable=Path("/env/bin/lifton"),
    )


def _git(root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *arguments],
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.strip()


def _clean_source_repo(tmp_path: Path) -> release.SourceSpec:
    root = tmp_path / "source"
    package = root / "lifton"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("VALUE = 1\n", encoding="utf-8")
    (root / ".gitignore").write_text(
        "__pycache__/\n*.py[cod]\nignored_probe.py\n",
        encoding="utf-8",
    )
    _git(root, "init", "--quiet")
    _git(root, "config", "user.name", "LiftOn Tests")
    _git(root, "config", "user.email", "lifton-tests@example.invalid")
    _git(root, "add", ".gitignore", "lifton/__init__.py")
    _git(root, "commit", "--quiet", "-m", "fixture")
    return release.SourceSpec(
        label="candidate",
        root=root,
        sha=_git(root, "rev-parse", "HEAD"),
        lifton_executable=Path(sys.executable).with_name("lifton"),
    )


def test_semantic_hash_ignores_record_and_attribute_order(tmp_path):
    first = tmp_path / "first.gff3"
    second = tmp_path / "second.gff3"
    first.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t12\t.\t+\t.\tID=g1;Name=one\n"
        "chr1\tLiftOn\tmRNA\t1\t12\t.\t+\t.\tParent=g1;ID=t1\n"
    )
    second.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tLiftOn\tgene\t1\t12\t.\t+\t.\tName=one;ID=g1\n"
    )

    left = release.gff3_fingerprints(first)
    right = release.gff3_fingerprints(second)

    assert left["byte_sha256"] != right["byte_sha256"]
    assert left["semantic_sha256"] == right["semantic_sha256"]
    assert left["semantic_algorithm"] == release.SEMANTIC_HASH_ALGORITHM
    assert left["feature_counts"] == {"gene": 1, "mRNA": 1}


def test_semantic_hash_is_reordered_multiset_and_duplicate_sensitive(tmp_path):
    rows = {
        "a": "chr1\tLiftOn\tgene\t1\t2\t.\t+\t.\tID=a;Name=A\n",
        "b": "chr1\tLiftOn\tgene\t3\t4\t.\t+\t.\tID=b;Name=B\n",
    }
    first = tmp_path / "first.gff3"
    reordered = tmp_path / "reordered.gff3"
    different_duplicates = tmp_path / "different-duplicates.gff3"
    first.write_text(
        "##gff-version 3\n" + rows["a"] + rows["b"] + rows["b"],
        encoding="utf-8",
    )
    reordered.write_text(
        "##gff-version 3\n" + rows["b"] + rows["a"] + rows["b"],
        encoding="utf-8",
    )
    different_duplicates.write_text(
        "##gff-version 3\n" + rows["a"] + rows["a"] + rows["b"],
        encoding="utf-8",
    )

    first_hash = release.gff3_fingerprints(first)
    reordered_hash = release.gff3_fingerprints(reordered)
    duplicates_hash = release.gff3_fingerprints(different_duplicates)

    assert first_hash["semantic_sha256"] == reordered_hash["semantic_sha256"]
    assert first_hash["semantic_sha256"] != duplicates_hash["semantic_sha256"]
    assert first_hash["feature_records"] == duplicates_hash["feature_records"] == 3


def test_semantic_multiset_accumulator_keeps_constant_row_state():
    accumulator = release._SemanticMultisetHash()

    for index in range(10_000):
        accumulator.add(f"row-{index}")

    assert accumulator.count == 10_000
    assert not hasattr(accumulator, "__dict__")
    assert isinstance(accumulator._sum, int)
    assert isinstance(accumulator._sum_squares, int)
    assert len(accumulator.hexdigest()) == 64


def test_stable_id_preservation_uses_declared_unique_same_type_ids(tmp_path):
    reference = tmp_path / "reference.gff3"
    output = tmp_path / "output.gff3"
    reference.write_text(
        "##gff-version 3\n"
        "chr1\tRef\tCDS\t1\t10\t.\t+\t0\tID=cds-A%3A1;Parent=tx%3A1\n"
        "chr1\tRef\tCDS\t21\t30\t.\t+\t0\tID=cds-A%3A1;Parent=tx%3A1\n"
        "chr1\tRef\texon\t1\t10\t.\t+\t.\tID=exon-A;Parent=tx1\n"
        "chr1\tRef\texon\t21\t30\t.\t+\t.\tParent=tx1\n",
        encoding="utf-8",
    )
    output.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tCDS\t101\t110\t.\t+\t0\tID=cds-A%3A1_1;Parent=tx%3A1_1\n"
        "chr1\tLiftOn\tCDS\t121\t130\t.\t+\t0\tID=cds-A%3A1_1;Parent=tx%3A1_1\n"
        # The exon ID under the wrong feature type must not count.
        "chr1\tLiftOn\tgene\t101\t130\t.\t+\t.\tID=exon-A\n",
        encoding="utf-8",
    )

    result = release.stable_id_preservation(reference, output)

    cds = result["by_type"]["CDS"]
    assert cds["n_reference_records"] == 2
    assert cds["n_reference_ids"] == 1
    assert cds["n_preserved_ids"] == 1
    assert cds["n_output_records"] == 2
    assert cds["n_output_records_with_id"] == 2
    assert cds["n_output_ids"] == 1
    assert cds["preservation_rate"] == 1.0
    exon = result["by_type"]["exon"]
    assert exon["n_reference_records"] == 2
    assert exon["n_reference_records_with_id"] == 1
    assert exon["n_reference_ids"] == 1
    assert exon["n_preserved_ids"] == 0
    assert exon["n_output_records"] == 0
    assert exon["n_output_records_with_id"] == 0
    assert exon["n_output_ids"] == 0
    assert exon["preservation_rate"] == 0.0


def test_stable_id_preservation_rejects_unrelated_numeric_suffix(tmp_path):
    reference = tmp_path / "reference.gff3"
    output = tmp_path / "output.gff3"
    reference.write_text(
        "##gff-version 3\n"
        "chr1\tRef\tCDS\t1\t10\t.\t+\t0\tID=cds-A;Parent=tx-A\n",
        encoding="utf-8",
    )
    output.write_text(
        "##gff-version 3\n"
        # A numeric ID suffix is not a copy marker unless Parent carries that
        # same suffix and maps back to the reference parent.
        "chr1\tLiftOn\tCDS\t101\t110\t.\t+\t0\t"
        "ID=cds-A_7;Parent=unrelated_7\n",
        encoding="utf-8",
    )

    cds = release.stable_id_preservation(reference, output)["by_type"]["CDS"]

    assert cds["n_output_ids"] == 1
    assert cds["n_preserved_ids"] == 0
    assert cds["preservation_rate"] == 0.0


def test_stable_id_preservation_marks_absent_or_idless_types_not_applicable(
        tmp_path):
    reference = tmp_path / "reference.gff3"
    output = tmp_path / "output.gff3"
    reference.write_text(
        "##gff-version 3\n"
        "chr1\tRef\texon\t1\t10\t.\t+\t.\tParent=tx1\n",
        encoding="utf-8",
    )
    output.write_text("##gff-version 3\n", encoding="utf-8")

    result = release.stable_id_preservation(reference, output)["by_type"]

    assert result["CDS"]["applicable"] is False
    assert result["CDS"]["reason"] == "reference_feature_type_absent"
    assert result["exon"]["applicable"] is False
    assert result["exon"]["reason"] == "no_declared_reference_ids"
    assert result["exon"]["n_output_records"] == 0
    assert result["exon"]["n_output_records_with_id"] == 0
    assert result["exon"]["preservation_rate"] is None


def test_subset_and_full_commands_pin_the_intended_protocol(tmp_path):
    source = _source(tmp_path)
    subset = release.build_lifton_argv(
        source, _inputs(tmp_path / "subset", "subset"),
        tmp_path / "subset.gff3", threads=8,
    )
    full = release.build_lifton_argv(
        source, _inputs(tmp_path / "full", "full"),
        tmp_path / "full.gff3", threads=8,
    )

    assert subset[subset.index("-t") + 1] == "1"
    assert "--locus-pipeline" not in subset
    assert "-T" not in subset and "-P" not in subset
    assert full[full.index("-t") + 1] == "8"
    assert "--locus-pipeline" in full
    assert "-T" in full and "-P" in full
    assert "--no-miniprot-rescue" in subset
    assert "--no-miniprot-rescue" in full


def test_e2e_input_resolution_uses_the_requested_dataset_registry(
        tmp_path, monkeypatch):
    registry = tmp_path / "custom-datasets.json"
    registry.write_text(json.dumps({
        "datasets": [{
            "id": "custom",
            "species": "Custom species",
            "reference_fa": "https://example/ref.fa",
            "target_fa": "https://example/target.fa",
            "reference_gff": "https://example/ref.gff3",
            "target_gff": "https://example/published.gff3",
        }],
    }))
    data_root = tmp_path / "data"
    dataset_dir = data_root / "custom"
    dataset_dir.mkdir(parents=True)
    for name in ("ref.fa", "target.fa", "ref.gff3", "published.gff3"):
        (dataset_dir / name).write_text(name + "\n")
    monkeypatch.setattr(release.run_benchmarks, "DEFAULT_DATA_DIR", data_root)

    inputs = release.resolve_panel_inputs(
        "e2e", "custom", dataset_registry=registry,
    )

    assert inputs.species == "Custom species"
    assert inputs.ref_fa == dataset_dir / "ref.fa"
    assert inputs.tgt_fa == dataset_dir / "target.fa"
    assert inputs.truth_gff is None
    assert set(release.input_fingerprints(inputs)) == {
        "ref_fa", "ref_gff", "tgt_fa",
    }


def test_e2e_arm_does_not_stage_optional_published_truth(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path / "inputs", "e2e")
    truth = tmp_path / "inputs" / "published.gff3"
    truth.write_text("unused truth\n", encoding="utf-8")
    inputs = release.dataclasses.replace(inputs, truth_gff=truth)
    dataset = release.run_benchmarks.Dataset(
        id="demo",
        species="Demo",
        reference_fa="https://example/ref.fa",
        target_fa="https://example/target.fa",
        reference_gff="https://example/ref.gff3",
        target_gff="https://example/published.gff3",
    )
    registry = release.run_benchmarks.Registry(
        datasets=[dataset],
        evaluation_flags=["-E"],
    )
    monkeypatch.setattr(
        release,
        "_e2e_dataset",
        lambda *_args: (registry, dataset),
    )
    monkeypatch.setattr(
        release,
        "e2e_flags",
        lambda *_args, **_kwargs: ["-t", "2"],
    )
    monkeypatch.setattr(
        release,
        "validate_e2e_biology",
        lambda _row: {"reference_features": 1},
    )
    monkeypatch.setattr(
        release,
        "_validation_document",
        lambda _path: {"is_valid": True, "n_errors": 0},
    )
    captured = {}

    def fake_run_dataset(selected, *, data_dir, results_dir, **_kwargs):
        captured["target_gff"] = selected.target_gff
        captured["data_dir"] = data_dir
        output = results_dir / selected.id / "lifton.gff3"
        output.parent.mkdir(parents=True)
        output.write_text(
            "##gff-version 3\n"
            "chr1\tLiftOn\tgene\t1\t2\t.\t+\t.\tID=g1\n",
            encoding="utf-8",
        )
        profile = {
            "exit_code": 0,
            "wall_clock_seconds": 1.0,
            "peak_rss_mb": 2.0,
        }
        return {
            "out_gff": str(output),
            "lift_profile": profile,
            "eval_profile": profile,
            "score_summary": {},
            "evaluation_summary": {},
        }

    monkeypatch.setattr(
        release.run_benchmarks, "run_dataset", fake_run_dataset,
    )
    version_dir = tmp_path / "arm"

    release._run_e2e_one(
        _source(tmp_path),
        inputs,
        version_dir,
        threads=2,
        mode="safe",
        dataset_registry=tmp_path / "datasets.json",
    )

    assert captured["target_gff"] is None
    assert not (
        captured["data_dir"] / "demo" / "published.gff3"
    ).exists()
    assert truth.read_text(encoding="utf-8") == "unused truth\n"


def test_full_input_resolution_uses_the_requested_benchmark_registry(
        tmp_path, monkeypatch):
    paths = {
        name: tmp_path / name
        for name in (
            "ref.fa", "ref.gff3", "target.fa",
            "liftoff.gff3", "miniprot.gff3",
            "transcripts.fa", "proteins.fa",
        )
    }
    for path in paths.values():
        path.write_text(path.name + "\n")
    registry = tmp_path / "custom-benchmarks.json"
    registry.write_text(json.dumps({
        "benchmarks": [{
            "id": "custom",
            "species": "Custom pair",
            "cross_species": True,
            "annotation_database": "Ensembl",
            "ref_genome": str(paths["ref.fa"]),
            "ref_gff": str(paths["ref.gff3"]),
            "tgt_genome": str(paths["target.fa"]),
        }],
    }))
    monkeypatch.setattr(
        release.devel_refresh,
        "_cached_inputs",
        lambda _benchmark: (
            paths["liftoff.gff3"], paths["miniprot.gff3"],
            paths["transcripts.fa"], paths["proteins.fa"],
        ),
    )

    inputs = release.resolve_panel_inputs(
        "full", "custom", benchmark_registry=registry,
    )

    assert inputs.species == "Custom pair"
    assert inputs.cross_species is True
    assert inputs.annotation_database == "Ensembl"
    assert inputs.ref_fa == paths["ref.fa"]
    assert inputs.tgt_fa == paths["target.fa"]


def test_paired_order_alternates_ab_ba(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path / "inputs", "subset")
    candidate = release.SourceSpec(
        "candidate", tmp_path, "a" * 40, Path("/env/bin/lifton"),
    )
    reference = release.SourceSpec(
        "reference", tmp_path, "b" * 40, Path("/env/bin/lifton"),
    )
    seen = []

    monkeypatch.setattr(release, "verify_source", lambda source: {"sha": source.sha})
    monkeypatch.setattr(
        release, "resolve_panel_inputs",
        lambda *_args, **_kwargs: inputs,
    )
    monkeypatch.setattr(release, "input_fingerprints", lambda _inputs: {})

    def fake_run(source, _inputs, version_dir, *, threads):
        seen.append(source.label)
        version_dir.mkdir(parents=True)
        output = version_dir / f"{source.label}.gff3"
        output.write_text(
            "##gff-version 3\n"
            f"chr1\tLiftOn\tgene\t1\t2\t.\t+\t.\tID={source.label}\n"
        )
        return output, {
            "profile": {
                "wall_clock_seconds": 2.0 if source.label == "candidate" else 4.0,
                "peak_rss_mb": 10.0 if source.label == "candidate" else 20.0,
            },
        }

    def fake_score(_inputs, _outputs, documents, _cell_dir, *, threads):
        for document in documents.values():
            document["summary"] = {"completeness_coding": 1.0}

    monkeypatch.setattr(release, "_run_one", fake_run)
    monkeypatch.setattr(release, "_score_pair", fake_score)

    first = release.run_paired_cell(
        panel="subset", benchmark="demo", repetition=1,
        candidate=candidate, reference=reference,
        cell_dir=tmp_path / "rep1", threads=8,
    )
    second = release.run_paired_cell(
        panel="subset", benchmark="demo", repetition=2,
        candidate=candidate, reference=reference,
        cell_dir=tmp_path / "rep2", threads=8,
    )

    assert seen == ["reference", "candidate", "candidate", "reference"]
    assert first["order"] == ["reference", "candidate"]
    assert second["order"] == ["candidate", "reference"]
    assert first["ratios"] == {"wall": 0.5, "peak_rss": 0.5}
    assert first["registries"] == {
        "benchmark": str(release.DEFAULT_BENCHMARK_REGISTRY.resolve()),
        "dataset": str(release.run_benchmarks.DEFAULT_REGISTRY.resolve()),
    }
    assert json.loads((tmp_path / "rep2" / "pair_result.json").read_text())[
        "repetition"
    ] == 2


def test_neutral_scoring_records_exact_transcript_evidence(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path / "inputs", "subset")
    outputs = {}
    documents = {}
    for label in ("candidate", "reference"):
        output = tmp_path / f"{label}.gff3"
        output.write_text("##gff-version 3\n", encoding="utf-8")
        outputs[label] = output
        documents[label] = {
            "profile": {
                "wall_clock_seconds": 1.0,
                "peak_rss_mb": 2.0,
            },
        }
    monkeypatch.setattr(
        release.evaluator,
        "build_reference",
        lambda *_args, **_kwargs: ({}, {}),
    )

    def fake_evaluate(tool, *_args, **kwargs):
        assert kwargs["ref_index"] == {}
        path = tmp_path / "cell" / "evaluation" / f"{tool}.transcripts.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "ref_mrna_id\trecovered\n"
            f"{tool}-tx\t1\n",
            encoding="utf-8",
        )
        return {
            "tool": tool,
            "benchmark": "demo",
            "transcripts_tsv": str(path),
        }

    monkeypatch.setattr(release.evaluator, "evaluate_tool", fake_evaluate)

    release._score_pair(
        inputs,
        outputs,
        documents,
        tmp_path / "cell",
        threads=2,
    )

    for label in ("candidate", "reference"):
        record = documents[label]["evaluation_artifacts"]["transcripts_tsv"]
        path = tmp_path / "cell" / "evaluation" / f"{label}.transcripts.tsv"
        assert record == {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": release.sha256_file(path),
        }
        assert documents[label]["summary"]["transcripts_tsv"] == str(
            path.resolve()
        )
        stable_ids = documents[label]["summary"]["stable_id_preservation"]
        assert stable_ids["method"] == release.STABLE_ID_METHOD
        assert all(
            row["applicable"] is False
            for row in stable_ids["by_type"].values()
        )


def test_neutral_scoring_never_creates_canonical_input_sidecars(
        tmp_path, monkeypatch):
    inputs = _inputs(tmp_path / "canonical", "subset")
    canonical = inputs.required_paths()
    before = {
        name: (path.read_bytes(), path.stat().st_mtime_ns)
        for name, path in canonical.items()
    }
    before_names = sorted(path.name for path in (tmp_path / "canonical").iterdir())
    observed_paths = []

    def fake_build_reference(ref_gff, ref_fa, **_kwargs):
        observed_paths.extend([Path(ref_gff), Path(ref_fa)])
        Path(str(ref_gff) + ".eval_db").write_text("local db\n")
        Path(str(ref_fa) + ".fai").write_text("local index\n")
        return {}, {}

    monkeypatch.setattr(
        release.evaluator, "build_reference", fake_build_reference,
    )

    def evaluate(tool, output, tgt_fa, reference, manifest, out_dir, *_args,
                 **_kwargs):
        del output, reference, manifest
        observed_paths.append(Path(tgt_fa))
        Path(str(tgt_fa) + ".fai").write_text("local index\n")
        transcript = Path(out_dir) / f"{tool}.transcripts.tsv"
        transcript.parent.mkdir(parents=True, exist_ok=True)
        transcript.write_text("ref_mrna_id\trecovered\n", encoding="utf-8")
        return {
            "tool": tool,
            "benchmark": "demo",
            "transcripts_tsv": str(transcript),
        }

    monkeypatch.setattr(release.evaluator, "evaluate_tool", evaluate)
    for repetition in (1, 2):
        cell = tmp_path / f"repetition-{repetition}"
        outputs = {}
        documents = {}
        for label in ("candidate", "reference"):
            output = cell / f"{label}.gff3"
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text("##gff-version 3\n", encoding="utf-8")
            outputs[label] = output
            documents[label] = {
                "profile": {
                    "wall_clock_seconds": 1.0,
                    "peak_rss_mb": 2.0,
                },
            }
        release._score_pair(
            inputs, outputs, documents, cell, threads=2,
        )

    assert sorted(
        path.name for path in (tmp_path / "canonical").iterdir()
    ) == before_names
    for name, path in canonical.items():
        assert (path.read_bytes(), path.stat().st_mtime_ns) == before[name]
        assert not Path(str(path) + ".fai").exists()
        assert not Path(str(path) + ".eval_db").exists()
    assert observed_paths
    assert all(
        "score-input" in path.parts
        for path in observed_paths
    )
    for repetition in (1, 2):
        score_root = tmp_path / f"repetition-{repetition}" / "score-input"
        assert (score_root / "reference-genome.fai").is_file()
        assert (score_root / "target-genome.fai").is_file()
        assert (score_root / "reference-annotation.eval_db").is_file()


def test_bootstrap_geomean_is_deterministic_and_ratio_based():
    first = release.bootstrap_geomean_ratio(
        [0.5, 1.0, 2.0], replicates=500, seed=17,
    )
    second = release.bootstrap_geomean_ratio(
        [0.5, 1.0, 2.0], replicates=500, seed=17,
    )

    assert first == second
    assert first["estimate"] == pytest.approx(1.0)
    assert first["low"] <= first["estimate"] <= first["high"]


def test_source_environment_pins_executable_and_import_paths(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", "/usr/bin")
    monkeypatch.setenv("PYTHONPATH", "/shared/python")
    source = release.SourceSpec(
        "candidate", tmp_path, "a" * 40, Path("/env/bin/lifton"),
    )

    environment = source.environment()

    assert environment["PATH"].split(":")[:2] == ["/env/bin", "/usr/bin"]
    assert environment["PYTHONPATH"].split(":")[:2] == [
        str(tmp_path), "/shared/python",
    ]
    assert environment["PYTHONDONTWRITEBYTECODE"] == "1"


def test_verify_source_is_repeatable_without_writing_bytecode(tmp_path):
    source = _clean_source_repo(tmp_path)

    first = release.verify_source(source)
    second = release.verify_source(source)

    assert first == second
    assert first["imported_package"] == str(
        (source.root / "lifton" / "__init__.py").resolve()
    )
    assert list(source.root.rglob("*.pyc")) == []
    assert _git(source.root, "status", "--porcelain") == ""


@pytest.mark.parametrize(
    ("relative_path", "category"),
    (
        ("lifton/untracked_probe.py", "untracked"),
        ("lifton/ignored_probe.py", "ignored"),
    ),
)
def test_verify_source_rejects_unversioned_importable_files(
        tmp_path, relative_path, category):
    source = _clean_source_repo(tmp_path)
    path = source.root / relative_path
    path.write_text("VALUE = 2\n", encoding="utf-8")

    with pytest.raises(
        RuntimeError,
        match=rf"unversioned import-affecting.*{category}:.*{path.name}",
    ):
        release.verify_source(source)


def test_e2e_flags_offer_fast_and_disk_backed_safe_modes(monkeypatch):
    registry = release.run_benchmarks.Registry(
        datasets=[],
        lifton_flags=[
            "--stream", "--inmemory-liftoff", "--locus-pipeline",
            "-t", "8", "--native", "-copies",
        ],
    )
    monkeypatch.setattr(
        release.run_benchmarks, "load_registry", lambda _path: registry,
    )

    fast = release.e2e_flags("fast", threads=4)
    safe = release.e2e_flags("safe", threads=2)
    stream = release.e2e_flags("stream", threads=3)
    inmemory_native = release.e2e_flags("inmemory-native", threads=1)

    assert fast[fast.index("-t") + 1] == "4"
    assert {"--stream", "--inmemory-liftoff", "--native"} <= set(fast)
    assert safe[safe.index("-t") + 1] == "2"
    assert not {"--stream", "--inmemory-liftoff", "--native"} & set(safe)
    assert "--locus-pipeline" in safe
    assert {"--stream"} == (
        {"--stream", "--inmemory-liftoff", "--native"} & set(stream)
    )
    assert {"--inmemory-liftoff", "--native"} == (
        {"--stream", "--inmemory-liftoff", "--native"} & set(inmemory_native)
    )
    with pytest.raises(ValueError, match="unknown end-to-end mode"):
        release.e2e_flags("turbo", threads=8)


def _e2e_row() -> dict:
    return {
        "biological_summary": {
            "schema_version": 1,
            "reference_features": 10,
            "mapped_features": 8,
            "lost_features": 2,
            "extra_copy_features": 1,
            "feature_completeness": 0.8,
            "emitted_transcript_records": 9,
            "mapped_transcripts_reported": 8,
            "evaluated_transcript_records": 9,
            "evaluated_coding_records": 7,
            "mean_protein_identity": 0.95,
        },
        "score_summary": {
            "format": "lifton_score_v1",
            "records": 9,
            "coding_records": 7,
            "noncoding_records": 2,
            "malformed_records": 0,
        },
        "evaluation_summary": {
            "format": "lifton_eval_v1",
            "records": 9,
            "coding_records": 7,
            "noncoding_records": 2,
            "malformed_records": 0,
        },
    }


def test_e2e_biology_validation_accepts_consistent_nonempty_schema():
    biological = release.validate_e2e_biology(_e2e_row())

    assert biological["feature_completeness"] == pytest.approx(0.8)
    assert biological["mean_protein_identity"] == pytest.approx(0.95)


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("biological_summary", "lost_features"), 3, "does not equal"),
        (("biological_summary", "mean_protein_identity"), float("nan"), "finite"),
        (("evaluation_summary", "records"), 0, "positive integer"),
        (("score_summary", "malformed_records"), 1, "not zero"),
    ],
)
def test_e2e_biology_validation_rejects_silent_bad_results(path, value, message):
    row = _e2e_row()
    row[path[0]][path[1]] = value

    with pytest.raises(RuntimeError, match=message):
        release.validate_e2e_biology(row)


def test_run_pair_cli_exposes_e2e_modes(tmp_path, monkeypatch):
    captured = {}
    benchmark_registry = tmp_path / "benchmarks.json"
    dataset_registry = tmp_path / "datasets.json"
    monkeypatch.setattr(
        release,
        "run_paired_cell",
        lambda **arguments: captured.update(arguments),
    )

    result = release.main([
        "run-pair",
        "--panel", "e2e",
        "--benchmark", "bee",
        "--repetition", "1",
        "--candidate-root", str(tmp_path / "candidate"),
        "--candidate-sha", "c" * 40,
        "--reference-root", str(tmp_path / "reference"),
        "--reference-sha", "r" * 40,
        "--lifton-executable", "/env/bin/lifton",
        "--cell-dir", str(tmp_path / "cell"),
        "--benchmark-registry", str(benchmark_registry),
        "--dataset-registry", str(dataset_registry),
        "--candidate-e2e-mode", "fast",
        "--reference-e2e-mode", "safe",
    ])

    assert result == 0
    assert captured["panel"] == "e2e"
    assert captured["candidate_e2e_mode"] == "fast"
    assert captured["reference_e2e_mode"] == "safe"
    assert captured["benchmark_registry"] == benchmark_registry.resolve()
    assert captured["dataset_registry"] == dataset_registry.resolve()
    assert captured["threads"] == 8


@pytest.mark.parametrize("threads", ("0", "-1"))
def test_run_pair_cli_rejects_nonpositive_threads(tmp_path, threads):
    arguments = [
        "run-pair",
        "--panel", "subset",
        "--benchmark", "demo",
        "--repetition", "1",
        "--candidate-root", str(tmp_path / "candidate"),
        "--candidate-sha", "c" * 40,
        "--reference-root", str(tmp_path / "reference"),
        "--reference-sha", "r" * 40,
        "--cell-dir", str(tmp_path / "cell"),
        "--threads", threads,
    ]

    with pytest.raises(SystemExit) as error:
        release.build_parser().parse_args(arguments)

    assert error.value.code == 2


@pytest.mark.parametrize("repetition", ("0", "6", "999"))
def test_run_pair_cli_rejects_out_of_range_repetition(tmp_path, repetition):
    arguments = [
        "run-pair",
        "--panel", "subset",
        "--benchmark", "demo",
        "--repetition", repetition,
        "--candidate-root", str(tmp_path / "candidate"),
        "--candidate-sha", "c" * 40,
        "--reference-root", str(tmp_path / "reference"),
        "--reference-sha", "r" * 40,
        "--cell-dir", str(tmp_path / "cell"),
    ]

    with pytest.raises(SystemExit) as error:
        release.build_parser().parse_args(arguments)

    assert error.value.code == 2


def test_run_paired_cell_rejects_direct_out_of_range_protocol_values(tmp_path):
    source = _source(tmp_path)

    with pytest.raises(ValueError, match="repetition must be between 1 and 5"):
        release.run_paired_cell(
            panel="subset",
            benchmark="demo",
            repetition=6,
            candidate=source,
            reference=source,
            cell_dir=tmp_path / "repetition",
        )
    with pytest.raises(ValueError, match="threads must be positive"):
        release.run_paired_cell(
            panel="subset",
            benchmark="demo",
            repetition=1,
            candidate=source,
            reference=source,
            cell_dir=tmp_path / "threads",
            threads=0,
        )


def test_release_arm_manifest_supports_pre_manifest_reference(tmp_path):
    version_dir = tmp_path / "reference"
    version_dir.mkdir()
    output = version_dir / "reference.gff3"
    output.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t2\t.\t+\t.\tID=g1\n"
    )
    native = version_dir / "lifton_output" / "run_manifest.json"
    fingerprints = release.gff3_fingerprints(output)

    path, document = release._write_release_arm_manifest(
        version_dir,
        source={"label": "reference", "sha": "a" * 40},
        protocol={"kind": "subset", "argv": ["lifton"]},
        profile={"wall_clock_seconds": 1.0, "peak_rss_mb": 2.0},
        output=output,
        fingerprints=fingerprints,
        validation={"is_valid": True, "n_errors": 0},
        native_manifests={"lift": native},
    )

    assert path == version_dir / "release_run_manifest.json"
    assert document["run"]["status"] == "success"
    assert document["artifacts"]["output_gff"]["byte_sha256"] == (
        fingerprints["byte_sha256"]
    )
    assert document["artifacts"]["native_manifests"]["lift"] == {
        "path": str(native.resolve()),
        "present": False,
    }
    assert json.loads(path.read_text())["source"]["sha"] == "a" * 40
    assert not path.with_suffix(".json.tmp").exists()
