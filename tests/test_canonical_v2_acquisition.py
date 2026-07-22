"""Security and readiness tests for canonical-v2 acquisition."""
from __future__ import annotations

import copy
import gzip
import json
import sys
import zipfile
from pathlib import Path

import pytest

from benchmarks import canonical_v2_acquisition as acquisition
from benchmarks import manifest_tools
from benchmarks import run_benchmarks


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = (
    REPO_ROOT / "benchmarks" / "manifests" / "canonical_v2_datasets.json"
)


@pytest.fixture(scope="module")
def manifest():
    return manifest_tools.validate_manifest(
        manifest_tools.load_json(MANIFEST_PATH)
    )


def _fasta() -> bytes:
    return b">chr22\n" + b"A" * 80 + b"\n"


def _gff() -> bytes:
    records = [b"##gff-version 3\n"]
    for index in range(100):
        records.extend([
            f"chr22\ttest\tgene\t1\t80\t.\t+\t.\tID=g{index}\n".encode(),
            (
                f"chr22\ttest\tmRNA\t1\t80\t.\t+\t.\t"
                f"ID=t{index};Parent=g{index}\n"
            ).encode(),
        ])
    return b"".join(records)


def _gtf() -> bytes:
    records = []
    for index in range(100):
        records.extend([
            (
                f'chr22\ttest\tgene\t1\t80\t.\t+\t.\t'
                f'gene_id "g{index}";\n'
            ).encode(),
            (
                f'chr22\ttest\ttranscript\t1\t80\t.\t+\t.\t'
                f'gene_id "g{index}"; transcript_id "t{index}";\n'
            ).encode(),
        ])
    return b"".join(records)


def _payload(role: str, filename: str) -> bytes:
    if "gtf" in filename:
        payload = _gtf()
    elif role in {"genome", "reference_genome", "target_genome"}:
        payload = _fasta()
    elif role in {"protein", "reference_proteins"}:
        payload = b">p1\nMAAA\n"
    else:
        payload = _gff()
    if filename.endswith(".gz"):
        return gzip.compress(payload, mtime=0)
    return payload


def _runner(manifest, calls):
    by_accession = {
        source["identity"]["accession"]: source
        for source in manifest["sources"]
        if source["transport"] == "ncbi_datasets"
    }

    def run(argv):
        argv = list(argv)
        calls.append(argv)
        if argv[0] == "curl":
            destination = Path(argv[argv.index("--output") + 1])
            url = argv[-1]
            filename = url.rsplit("/", 1)[-1]
            source = next(
                source
                for source in manifest["sources"]
                if any(
                    item["locator"].get("filename") == filename
                    for item in source["files"]
                )
            )
            role = next(
                item["role"] for item in source["files"]
                if item["locator"]["filename"] == filename
            )
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(_payload(role, filename))
            return 0
        accession = argv[argv.index("accession") + 1]
        source = by_accession[accession]
        destination = Path(argv[argv.index("--filename") + 1])
        destination.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(destination, "w") as handle:
            for item in source["files"]:
                member = item["locator"]["archive_member"]
                handle.writestr(
                    member,
                    _payload(item["role"], Path(member).name),
                )
        return 0

    return run


def _existing_registry(tmp_path):
    entries = []
    for benchmark_id in ("arabidopsis", "bee", "human_mane"):
        root = tmp_path / "existing" / benchmark_id
        root.mkdir(parents=True)
        paths = {
            "ref_genome": root / "reference.fa",
            "ref_gff": root / "reference.gff3",
            "ref_proteins": root / "reference.faa",
            "tgt_genome": root / "target.fa",
        }
        paths["ref_genome"].write_bytes(_fasta())
        paths["tgt_genome"].write_bytes(_fasta())
        paths["ref_gff"].write_bytes(_gff())
        paths["ref_proteins"].write_bytes(b">p1\nMAAA\n")
        entries.append({
            "id": benchmark_id,
            "species": benchmark_id,
            "cross_species": False,
            "ref_genome": str(paths["ref_genome"]),
            "ref_gff": str(paths["ref_gff"]),
            "ref_proteins": (
                None if benchmark_id == "human_mane"
                else str(paths["ref_proteins"])
            ),
            "tgt_genome": str(paths["tgt_genome"]),
            "ref_chrom": "chr22",
            "tgt_chrom": "chr22",
            "annotation_database": "RefSeq",
        })
    registry = tmp_path / "benchmarks.json"
    registry.write_text(json.dumps({
        "tools": {
            "samtools_bin": "samtools",
            "minimap2_bin": "minimap2",
        },
        "benchmarks": entries,
    }))
    return registry


def _dataset_registry(tmp_path):
    path = tmp_path / "datasets.json"
    path.write_text(json.dumps({
        "datasets": [],
        "lifton_flags": [],
        "evaluation_flags": [],
    }))
    return path


def _synthetic_builder(
    scenario,
    source_fasta,
    source_gff,
    output_dir,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / "target.fa"
    truth = output_dir / "target.truth.gff3"
    ortholog = output_dir / "ortholog_map.json"
    target.write_bytes(_fasta())
    truth.write_bytes(_gff())
    ortholog.write_text(json.dumps({
        "mappings": [
            {"source_id": "g1", "truth_ids": ["g1"], "feature_type": "gene"},
        ],
    }))
    outputs = {
        "target_fasta": target,
        "truth_gff": truth,
        "ortholog_map": ortholog,
    }
    transform = {
        "schema_version": 1,
        "kind": "lifton-synthetic-target-truth",
        "source": {
            "seqid": scenario["design"]["chromosome"],
            "fasta": {"sha256": manifest_tools.sha256_file(source_fasta)},
            "gff": {"sha256": manifest_tools.sha256_file(source_gff)},
        },
        "transform": manifest_tools.expected_synthetic_transform(scenario),
        "outputs": {
            role: {
                "name": path.name,
                "sha256": manifest_tools.sha256_file(path),
            }
            for role, path in outputs.items()
        },
    }
    manifest_path = output_dir / "transform.manifest.json"
    manifest_path.write_text(json.dumps(transform))
    return {**outputs, "manifest": manifest_path}


def _ortholog_registry(tmp_path, manifest, resolved_scenarios):
    records = {}
    input_roles = {
        "source_annotation": "reference_annotation",
        "source_genome": "reference_genome",
        "target_annotation": "target_truth",
        "target_genome": "target_genome",
    }
    for scenario_id in acquisition.BIOLOGICAL_IDS:
        inputs = {
            provenance_role: Path(
                resolved_scenarios[scenario_id][scenario_role]
            )
            for provenance_role, scenario_role in input_roles.items()
        }
        mappings = [
            {
                "source_id": f"g{index}",
                "truth_ids": [f"g{index}"],
                "feature_type": "gene",
                "status": "retained",
            }
            for index in range(100)
        ] + [
            {
                "source_id": f"t{index}",
                "truth_ids": [f"t{index}"],
                "feature_type": "transcript",
                "status": "retained",
            }
            for index in range(100)
        ]
        mapping = tmp_path / f"{scenario_id}.orthologs.json"
        mapping.write_text(json.dumps({
            "schema_version": 1,
            "method": "protein-rbh-ortholog-scope-v1",
            "metadata": {
                "scope": "full",
                "parameters": {},
                "provenance": {
                    "inputs": {
                        role: {
                            "path": str(input_path.resolve()),
                            "size": input_path.stat().st_size,
                            "sha256": manifest_tools.sha256_file(input_path),
                        }
                        for role, input_path in inputs.items()
                    },
                    "tools": {
                        "gffread": {"version": "fixture"},
                        "mmseqs": {"version": "fixture"},
                    },
                    "commands": [{
                        "argv": ["fixture-rbh"],
                        "shell": False,
                    }],
                },
                "counts": {
                    "source_genes": 100,
                    "source_transcripts": 100,
                    "target_genes": 100,
                    "target_transcripts": 100,
                    "gene_groups_retained": 100,
                    "transcript_groups_retained": 100,
                    "gene_groups_unscored": 0,
                    "transcript_groups_unscored": 0,
                },
            },
            "mappings": mappings,
        }))
        records[scenario_id] = {
            "path": mapping.name,
            "sha256": manifest_tools.sha256_file(mapping),
            "id_policy": "ortholog-map",
        }
    path = tmp_path / "ortholog-registry.json"
    path.write_text(json.dumps({
        "schema_version": 1,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "mappings": records,
    }))
    return path


def _prepared_subset_manifest(work_root: Path, scenario_id: str) -> None:
    subset = work_root / scenario_id / "subset"
    subset.mkdir(parents=True)
    source_gff = subset / "ref.chrom.gff3"
    target_fasta = subset / "tgt.chrom.fa"
    source_gff.write_bytes(_gff())
    target_fasta.write_bytes(_fasta())
    document = {
        "schema_version": 2,
        "id": scenario_id,
        "paths": {
            "ref_gff": str(source_gff.resolve()),
            "tgt_fa": str(target_fasta.resolve()),
        },
        "outputs": {
            "ref_gff": acquisition._file_record(source_gff),
            "tgt_fa": acquisition._file_record(target_fasta),
        },
    }
    (subset / "subset.manifest.json").write_text(json.dumps(document))


def test_dry_run_is_non_mutating_and_names_first_observation_policy(
    manifest,
    tmp_path,
):
    cache_root = tmp_path / "cache"
    result = acquisition.acquisition_dry_run(manifest, cache_root)
    assert result["execution"]["shell"] is False
    assert result["execution"]["requires_identity_pin_acknowledgement"] is True
    assert result["execution"]["identity_pinned_roles"]
    commands = [
        command
        for step in result["steps"]
        for command in step["commands"]
    ]
    assert all(
        command[-2 if command[0] == "curl" else -1].endswith(".part")
        for command in commands
    )
    assert all(
        "--continue-at" in command
        for command in commands
        if command[0] == "curl"
    )
    assert not cache_root.exists()


def test_lock_path_is_contained_and_dry_run_reports_custom_path(
    manifest,
    tmp_path,
):
    cache_root = tmp_path / "cache"
    custom_lock = cache_root / "locks" / "canonical.lock.json"

    result = acquisition.acquisition_dry_run(
        manifest,
        cache_root,
        lock_path=custom_lock,
    )

    assert result["execution"]["atomic_lock"] == str(custom_lock.resolve())
    assert not cache_root.exists()
    outside = tmp_path / "outside.lock.json"
    with pytest.raises(acquisition.AcquisitionError, match="child of cache root"):
        acquisition.acquisition_dry_run(
            manifest,
            cache_root,
            lock_path=outside,
        )
    with pytest.raises(acquisition.AcquisitionError, match="child of cache root"):
        acquisition.acquire(
            manifest,
            cache_root=cache_root,
            lock_path=outside,
            accept_identity_pinned_bytes=True,
            runner=lambda _argv: pytest.fail("runner must not execute"),
        )
    assert not outside.exists()


def test_acquisition_requires_acknowledgement_locks_and_resumes(
    manifest,
    tmp_path,
):
    calls = []
    cache_root = tmp_path / "cache"
    with pytest.raises(
        acquisition.AcquisitionError,
        match="accept-identity-pinned-bytes",
    ):
        acquisition.acquire(
            manifest,
            cache_root=cache_root,
            runner=_runner(manifest, calls),
        )
    assert calls == []

    lock, result = acquisition.acquire(
        manifest,
        cache_root=cache_root,
        accept_identity_pinned_bytes=True,
        runner=_runner(manifest, calls),
    )
    assert result["verified"] is True
    assert result["resumed_from_verified_lock"] is False
    assert set(lock["sources"]) == {
        source["id"] for source in manifest["sources"]
        if source["transport"] != "existing_registry"
    }
    assert all(Path(record["path"]).parts[0] == "sha256"
               for source in lock["sources"].values()
               for record in source["files"].values())

    def must_not_run(_argv):
        raise AssertionError("verified acquisition must resume without commands")

    resumed_lock, resumed = acquisition.acquire(
        manifest,
        cache_root=cache_root,
        runner=must_not_run,
    )
    assert resumed_lock == lock
    assert resumed["resumed_from_verified_lock"] is True


def test_zip_validation_rejects_traversal_and_links(tmp_path):
    traversal = tmp_path / "traversal.zip"
    with zipfile.ZipFile(traversal, "w") as handle:
        handle.writestr("../escape", b"payload")
    with pytest.raises(acquisition.AcquisitionError, match="unsafe archive"):
        acquisition.validate_zip_archive(traversal)

    link = tmp_path / "link.zip"
    info = zipfile.ZipInfo("safe/link")
    info.create_system = 3
    info.external_attr = (stat_mode := 0o120777) << 16
    assert stat_mode
    with zipfile.ZipFile(link, "w") as handle:
        handle.writestr(info, b"target")
    with pytest.raises(acquisition.AcquisitionError, match="links"):
        acquisition.validate_zip_archive(link)


def test_declared_checksum_mismatch_fails_closed(manifest, tmp_path):
    changed = copy.deepcopy(manifest)
    item = changed["sources"][0]["files"][0]
    item["pin_state"] = "sha256_pinned"
    item["expected_sha256"] = "0" * 64
    changed = manifest_tools.validate_manifest(changed)
    with pytest.raises(acquisition.AcquisitionError, match="SHA-256 mismatch"):
        acquisition.acquire(
            changed,
            cache_root=tmp_path / "cache",
            accept_identity_pinned_bytes=True,
            runner=_runner(changed, []),
        )
    assert not (tmp_path / "cache" / acquisition.LOCK_NAME).exists()


def test_gzip_materialization_fingerprints_bytes_and_builds_fai(tmp_path):
    source = tmp_path / "genome.fa.gz"
    source.write_bytes(gzip.compress(_fasta(), mtime=0))
    runtime, record = acquisition._materialize_asset(
        source,
        tmp_path / "runtime",
        fasta=True,
    )
    assert runtime.read_bytes() == _fasta()
    assert record["transform"] == "gzip-decompress"
    assert record["source"]["sha256"] != record["runtime"]["sha256"]
    assert Path(record["fasta_index"]["path"]).read_text() == (
        "chr22\t80\t7\t80\t81\n"
    )


def test_synthetic_cache_binds_design_transform_sources_and_outputs(
    manifest,
    tmp_path,
):
    scenario = copy.deepcopy(manifest["scenarios"][8])
    source_fasta = tmp_path / "source.fa"
    source_gff = tmp_path / "source.gff3"
    source_fasta.write_bytes(_fasta())
    source_gff.write_bytes(_gff())
    output = tmp_path / "synthetic"
    _synthetic_builder(scenario, source_fasta, source_gff, output)

    verified = acquisition._verified_synthetic_outputs(
        output, scenario, source_fasta, source_gff,
    )
    assert verified is not None
    _paths, provenance = verified
    assert provenance["design_sha256"] == manifest_tools.canonical_sha256(
        scenario["design"]
    )
    assert set(provenance["outputs"]["target_fasta"]) == {
        "path", "size", "sha256",
    }

    changed = copy.deepcopy(scenario)
    changed["design"]["operations"][0]["cut_after"][0] += 1
    changed = {
        **changed,
        "design": dict(changed["design"]),
    }
    assert acquisition._verified_synthetic_outputs(
        output, changed, source_fasta, source_gff,
    ) is None

    original_source_gff = source_gff.read_bytes()
    source_gff.write_bytes(original_source_gff + b"# changed source\n")
    assert acquisition._verified_synthetic_outputs(
        output, scenario, source_fasta, source_gff,
    ) is None
    source_gff.write_bytes(original_source_gff)

    transform_path = output / "transform.manifest.json"
    transform_document = json.loads(transform_path.read_text(encoding="utf-8"))
    transform_document["transform"]["cuts_after_source_coordinate"][0] += 1
    transform_path.write_text(json.dumps(transform_document), encoding="utf-8")
    assert acquisition._verified_synthetic_outputs(
        output, scenario, source_fasta, source_gff,
    ) is None
    _synthetic_builder(scenario, source_fasta, source_gff, output)

    target = output / "target.fa"
    target.write_bytes(target.read_bytes() + b"A\n")
    assert acquisition._verified_synthetic_outputs(
        output, scenario, source_fasta, source_gff,
    ) is None


def test_chr22_scope_resolves_exact_ncbi_accession(tmp_path):
    fasta = tmp_path / "reference.fa"
    annotation = tmp_path / "reference.gff3"
    fasta.write_text(">NC_000022.11 chromosome 22\n" + "A" * 80 + "\n")
    annotation.write_text(
        "##gff-version 3\n"
        "NC_000022.11\tRefSeq\tregion\t1\t80\t.\t+\t."
        "\tID=NC_000022.11:1..80;Name=22;chromosome=22\n"
        "NC_000022.11\tRefSeq\tgene\t1\t80\t.\t+\t.\tID=g1\n"
    )
    assert acquisition._reference_subset_seqid({
        "reference_genome": fasta,
        "reference_annotation": annotation,
    }) == "NC_000022.11"


def test_materialization_blocks_missing_curated_ortholog_maps(
    manifest,
    tmp_path,
):
    cache_root = tmp_path / "cache"
    lock, _result = acquisition.acquire(
        manifest,
        cache_root=cache_root,
        accept_identity_pinned_bytes=True,
        runner=_runner(manifest, []),
    )
    report = acquisition.materialize(
        manifest,
        lock,
        cache_root=cache_root,
        benchmark_registry=_existing_registry(tmp_path),
        dataset_registry=_dataset_registry(tmp_path),
        data_root=tmp_path / "data",
        work_root=tmp_path / "work",
        synthetic_builder=_synthetic_builder,
    )
    missing = [
        blocker for blocker in report["blockers"]
        if blocker["code"] == "missing_ortholog_map"
    ]
    assert [row["scenario_id"] for row in missing] == sorted(
        acquisition.BIOLOGICAL_IDS
    )
    assert report["campaign_ready"] is False
    assert report["registries_exported"] is False
    assert report["synthetic_provenance"]
    assert not (
        cache_root / "runtime" / "registries" / "benchmarks.json"
    ).exists()


def test_failed_rematerialization_invalidates_previous_ready_preflight(
        manifest, tmp_path):
    cache_root = tmp_path / "cache"
    cache_root.mkdir()
    preflight = cache_root / acquisition.PREFLIGHT_NAME
    preflight.write_text(json.dumps({
        "schema_version": 2,
        "campaign_ready": True,
        "registries_exported": True,
        "blockers": [],
        "remaining_actions": [],
    }))

    with pytest.raises(acquisition.AcquisitionError, match="benchmark registry"):
        acquisition.materialize(
            manifest,
            {},
            cache_root=cache_root,
            benchmark_registry=tmp_path / "missing-benchmarks.json",
            dataset_registry=tmp_path / "missing-datasets.json",
        )

    invalidated = json.loads(preflight.read_text())
    assert invalidated["campaign_ready"] is False
    assert invalidated["registries_exported"] is False
    assert invalidated["blockers"][0]["code"] == "materialization_in_progress"


def test_complete_maps_export_only_preparation_overlays_and_fail_closed(
    manifest,
    tmp_path,
):
    cache_root = tmp_path / "cache"
    lock, _result = acquisition.acquire(
        manifest,
        cache_root=cache_root,
        accept_identity_pinned_bytes=True,
        runner=_runner(manifest, []),
    )
    benchmark_registry = _existing_registry(tmp_path)
    dataset_registry = _dataset_registry(tmp_path)
    initial = acquisition.materialize(
        manifest,
        lock,
        cache_root=cache_root,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
        data_root=tmp_path / "data",
        work_root=tmp_path / "work",
        synthetic_builder=_synthetic_builder,
    )
    ortholog_registry = _ortholog_registry(
        tmp_path,
        manifest,
        initial["resolved_scenarios"],
    )
    report = acquisition.materialize(
        manifest,
        lock,
        cache_root=cache_root,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
        ortholog_registry=ortholog_registry,
        data_root=tmp_path / "data",
        work_root=tmp_path / "work",
        synthetic_builder=_synthetic_builder,
    )
    assert report["registries_exported"] is False
    assert report["preparation_registries_exported"] is True
    assert report["campaign_ready"] is False
    assert len(report["remaining_actions"]) == 10
    assert all("--registry" in argv for argv in report["remaining_actions"])
    assert {blocker["code"] for blocker in report["blockers"]} == {
        "subset_inputs_not_prepared",
    }
    assert report["schema_version"] == 2
    assert report["acquisition_lock_sha256"] == (
        manifest_tools.canonical_sha256(lock)
    )
    benchmark_overlay = Path(
        report["registries"]["benchmark_preparation"]["path"]
    )
    assert benchmark_overlay.name == "benchmarks.preparation.json"
    benchmark_document = json.loads(benchmark_overlay.read_text())
    generated = {
        row["id"]: row for row in benchmark_document["benchmarks"]
        if row["id"].startswith("v2_")
    }
    assert set(generated) == {
        *acquisition.BIOLOGICAL_IDS,
        "v2_synth_chr22_fragmented",
        "v2_synth_chr22_sv",
    }
    assert all(
        row["target_truth"]["id_policy"] == "ortholog-map"
        for row in generated.values()
    )
    assert generated["v2_deep_zebrafish_xenopus"]["tgt_chrom"] == "WHOLE"
    assert generated["v2_deep_tomato_rice"]["tgt_chrom"] == "WHOLE"
    assert generated["v2_synth_chr22_fragmented"]["tgt_chrom"] == "WHOLE"
    assert generated["v2_truth_human_grch38_chm13"][
        "full_input_mode"
    ] == "raw"

    dataset_overlay = Path(
        report["registries"]["dataset_preparation"]["path"]
    )
    registry = run_benchmarks.load_registry(dataset_overlay)
    assert {dataset.id for dataset in registry.datasets} == {
        "v2_dialect_flybase_dmel_dere",
        "v2_dialect_wormbase_celegans_cbriggsae",
        "v2_truth_rat_mouse_e2e",
    }
    for dataset in registry.datasets:
        assert dataset.cross_species is True
        assert dataset.truth_gff
        assert dataset.ortholog_map
        assert dataset.truth_id_policy == "ortholog-map"

    registry_document = json.loads(ortholog_registry.read_text())
    scenario_id = acquisition.BIOLOGICAL_IDS[0]
    mapping_record = registry_document["mappings"][scenario_id]
    mapping_path = ortholog_registry.parent / mapping_record["path"]
    mapping_document = json.loads(mapping_path.read_text())
    mapping_document["mappings"] = mapping_document["mappings"][:2]
    mapping_path.write_text(json.dumps(mapping_document))
    mapping_record["sha256"] = manifest_tools.sha256_file(mapping_path)
    ortholog_registry.write_text(json.dumps(registry_document))
    with pytest.raises(
        acquisition.AcquisitionError,
        match="omits source gene ID",
    ):
        acquisition.materialize(
            manifest,
            lock,
            cache_root=cache_root,
            benchmark_registry=benchmark_registry,
            dataset_registry=dataset_registry,
            ortholog_registry=ortholog_registry,
            data_root=tmp_path / "data",
            work_root=tmp_path / "work",
            synthetic_builder=_synthetic_builder,
        )


def test_complete_verified_inputs_export_campaign_ready_overlays(
    manifest,
    tmp_path,
    monkeypatch,
):
    cache_root = tmp_path / "cache"
    work_root = tmp_path / "work"
    lock, _result = acquisition.acquire(
        manifest,
        cache_root=cache_root,
        accept_identity_pinned_bytes=True,
        runner=_runner(manifest, []),
    )
    benchmark_registry = _existing_registry(tmp_path)
    dataset_registry = _dataset_registry(tmp_path)
    initial = acquisition.materialize(
        manifest,
        lock,
        cache_root=cache_root,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
        data_root=tmp_path / "data",
        work_root=work_root,
        synthetic_builder=_synthetic_builder,
    )
    ortholog_registry = _ortholog_registry(
        tmp_path,
        manifest,
        initial["resolved_scenarios"],
    )
    for scenario_id in acquisition.BIOLOGICAL_IDS:
        _prepared_subset_manifest(work_root, scenario_id)
    monkeypatch.setattr(
        acquisition,
        "_subset_preflight",
        lambda *_args, **_kwargs: ([], []),
    )

    report = acquisition.materialize(
        manifest,
        lock,
        cache_root=cache_root,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
        ortholog_registry=ortholog_registry,
        data_root=tmp_path / "data",
        work_root=work_root,
        synthetic_builder=_synthetic_builder,
    )

    assert report["schema_version"] == 2
    assert report["campaign_ready"] is True
    assert report["registries_exported"] is True
    assert report["blockers"] == []
    assert report["remaining_actions"] == []
    assert set(report["panel_truth_provenance"]) == set(
        acquisition.BIOLOGICAL_IDS
    )
    benchmark_overlay = Path(report["registries"]["benchmark"]["path"])
    rows = {
        row["id"]: row
        for row in json.loads(benchmark_overlay.read_text())["benchmarks"]
    }
    for scenario_id in acquisition.BIOLOGICAL_IDS:
        truth = rows[scenario_id]["target_truth_by_panel"]
        assert set(truth) == {"subset", "full"}
        assert Path(truth["subset"]["gff"]).is_file()
        assert Path(truth["subset"]["ortholog_map"]).is_file()


def test_bootstrap_resumes_actions_and_stops_only_at_ready(
    tmp_path,
    monkeypatch,
):
    from benchmarks.compare import ortholog_scope

    cache_root = tmp_path / "cache"
    cache_root.mkdir()
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text("{}")
    lock_path = cache_root / acquisition.LOCK_NAME
    lock_path.write_text("{}")
    preflight = cache_root / acquisition.PREFLIGHT_NAME
    benchmark_overlay = cache_root / "runtime" / "registries" / "benchmarks.json"
    dataset_overlay = cache_root / "runtime" / "registries" / "datasets.json"
    benchmark_overlay.parent.mkdir(parents=True)
    benchmark_overlay.write_text("{}")
    dataset_overlay.write_text("{}")
    reports = [
        {
            "campaign_ready": False,
            "blockers": [{"code": "missing_ortholog_map"}],
            "remaining_actions": [[
                sys.executable,
                "-B",
                "-m",
                "benchmarks.compare.ortholog_scope",
                "build",
            ]],
        },
        {
            "campaign_ready": False,
            "blockers": [{"code": "subset_inputs_not_prepared"}],
            "remaining_actions": [[
                sys.executable,
                "-B",
                "-m",
                "benchmarks.compare.build_inputs",
                "demo",
            ]],
        },
        {
            "campaign_ready": True,
            "blockers": [],
            "remaining_actions": [],
            "registries": {
                "benchmark": acquisition._file_record(benchmark_overlay),
                "dataset": acquisition._file_record(dataset_overlay),
            },
        },
    ]
    rounds = []

    monkeypatch.setattr(acquisition, "_safe_cache_root", lambda path: Path(path))
    monkeypatch.setattr(
        acquisition,
        "acquire",
        lambda *_args, **_kwargs: ({}, {"verified": True}),
    )

    def fake_materialize(*_args, **_kwargs):
        report = reports.pop(0)
        preflight.write_text(json.dumps(report))
        return report

    monkeypatch.setattr(acquisition, "materialize", fake_materialize)
    state_calls = []

    def fake_ortholog_state(*_args, ortholog_root, final_registry, **_kwargs):
        state_calls.append(True)
        if len(state_calls) == 3:
            Path(final_registry).parent.mkdir(parents=True, exist_ok=True)
            Path(final_registry).write_text("{}")
            return {
                "valid_scenarios": ["demo"],
                "missing_scenarios": [],
                "registry_path": Path(final_registry),
            }
        return {
            "valid_scenarios": [],
            "missing_scenarios": ["demo"],
            "registry_path": None,
        }

    monkeypatch.setattr(
        acquisition, "_bootstrap_ortholog_state", fake_ortholog_state,
    )
    monkeypatch.setattr(
        acquisition,
        "_resume_safe_actions",
        lambda actions, _state: actions,
    )
    monkeypatch.setattr(
        acquisition,
        "_run_bootstrap_actions",
        lambda actions, **_kwargs: rounds.append(actions) or [{
            "index": 1,
            "returncode": 0,
        }],
    )

    def fake_finalize(_manifest, _root, output):
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text("{}")

    monkeypatch.setattr(
        ortholog_scope, "finalize_mapping_registry", fake_finalize,
    )

    result = acquisition.bootstrap(
        {},
        manifest_path=manifest_path,
        cache_root=cache_root,
        lock_path=lock_path,
        benchmark_registry=tmp_path / "base-benchmarks.json",
        dataset_registry=tmp_path / "base-datasets.json",
        threads=8,
        max_active=2,
        accept_identity_pinned_bytes=True,
    )

    assert result["campaign_ready"] is True
    assert len(rounds) == 2
    assert len(result["action_rounds"]) == 2
    assert result["preflight"] == acquisition._file_record(preflight)


def test_bootstrap_finalizes_completed_bundles_before_materializing(
        tmp_path, monkeypatch):
    from benchmarks.compare import ortholog_scope

    scenario_ids = ("bio-one", "bio-two")
    manifest = {
        "scenarios": [
            {"id": scenario_id, "kind": "biological"}
            for scenario_id in scenario_ids
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    cache_root = tmp_path / "cache"
    (cache_root / acquisition.LOCK_NAME).parent.mkdir(parents=True)
    (cache_root / acquisition.LOCK_NAME).write_text("{}")
    ortholog_root = cache_root / "runtime" / "ortholog_scopes"
    for scenario_id in scenario_ids:
        bundle = ortholog_root / scenario_id
        bundle.mkdir(parents=True)
        (bundle / "ortholog_map.json").write_text(scenario_id)

    def fake_validate(bundle):
        mapping_path = Path(bundle) / "ortholog_map.json"
        return {
            "mapping_path": mapping_path,
            "mapping_sha256": manifest_tools.sha256_file(mapping_path),
        }

    def fake_finalize(_manifest_path, mapping_root, output):
        document = {
            "schema_version": acquisition.ORTHOLOG_SCHEMA_VERSION,
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "mappings": {
                scenario_id: {
                    "id_policy": "ortholog-map",
                    "path": f"{scenario_id}/ortholog_map.json",
                    "sha256": manifest_tools.sha256_file(
                        Path(mapping_root) / scenario_id / "ortholog_map.json"
                    ),
                }
                for scenario_id in scenario_ids
            },
        }
        Path(output).write_text(json.dumps(document))
        return document

    monkeypatch.setattr(ortholog_scope, "validate_scope_bundle", fake_validate)
    monkeypatch.setattr(
        ortholog_scope, "finalize_mapping_registry", fake_finalize,
    )
    monkeypatch.setattr(acquisition, "_safe_cache_root", lambda path: Path(path))
    monkeypatch.setattr(
        acquisition,
        "acquire",
        lambda *_args, **_kwargs: ({}, {"verified": True}),
    )
    preflight = cache_root / acquisition.PREFLIGHT_NAME
    overlay_root = cache_root / "runtime" / "registries"
    overlay_root.mkdir(parents=True)
    benchmark_overlay = overlay_root / "benchmarks.json"
    dataset_overlay = overlay_root / "datasets.json"
    benchmark_overlay.write_text("{}")
    dataset_overlay.write_text("{}")

    def fake_materialize(*_args, ortholog_registry, **_kwargs):
        assert Path(ortholog_registry).name == "ortholog_registry.json"
        report = {
            "campaign_ready": True,
            "blockers": [],
            "remaining_actions": [],
            "registries": {
                "benchmark": acquisition._file_record(benchmark_overlay),
                "dataset": acquisition._file_record(dataset_overlay),
            },
        }
        preflight.write_text(json.dumps(report))
        return report

    monkeypatch.setattr(acquisition, "materialize", fake_materialize)

    result = acquisition.bootstrap(
        manifest,
        manifest_path=manifest_path,
        cache_root=cache_root,
        lock_path=cache_root / acquisition.LOCK_NAME,
        benchmark_registry=tmp_path / "benchmarks.json",
        dataset_registry=tmp_path / "datasets.json",
        threads=1,
        max_active=1,
        accept_identity_pinned_bytes=True,
    )

    assert result["campaign_ready"] is True
    assert result["action_rounds"] == []
    assert result["ortholog_bootstrap"]["valid_scenarios"] == list(
        scenario_ids
    )
    assert (ortholog_root / "ortholog_registry.json").is_file()


def test_bootstrap_reuses_successes_after_partial_action_failure(
        tmp_path, monkeypatch):
    from benchmarks.compare import ortholog_scope

    scenario_ids = ("bio-one", "bio-two")
    manifest = {
        "scenarios": [
            {"id": scenario_id, "kind": "biological"}
            for scenario_id in scenario_ids
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    cache_root = tmp_path / "cache"
    (cache_root / acquisition.LOCK_NAME).parent.mkdir(parents=True)
    (cache_root / acquisition.LOCK_NAME).write_text("{}")
    ortholog_root = cache_root / "runtime" / "ortholog_scopes"
    preflight = cache_root / acquisition.PREFLIGHT_NAME
    overlay_root = cache_root / "runtime" / "registries"
    overlay_root.mkdir(parents=True)
    benchmark_overlay = overlay_root / "benchmarks.json"
    dataset_overlay = overlay_root / "datasets.json"
    benchmark_overlay.write_text("{}")
    dataset_overlay.write_text("{}")

    def publish_bundle(scenario_id):
        bundle = ortholog_root / scenario_id
        bundle.mkdir(parents=True)
        (bundle / "ortholog_map.json").write_text(scenario_id)

    def fake_validate(bundle):
        mapping_path = Path(bundle) / "ortholog_map.json"
        return {
            "mapping_path": mapping_path,
            "mapping_sha256": manifest_tools.sha256_file(mapping_path),
        }

    def registry_document(mapping_root, selected):
        return {
            "schema_version": acquisition.ORTHOLOG_SCHEMA_VERSION,
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "mappings": {
                scenario_id: {
                    "id_policy": "ortholog-map",
                    "path": f"{scenario_id}/ortholog_map.json",
                    "sha256": manifest_tools.sha256_file(
                        Path(mapping_root) / scenario_id / "ortholog_map.json"
                    ),
                }
                for scenario_id in selected
            },
        }

    def fake_finalize(_manifest_path, mapping_root, output):
        document = registry_document(mapping_root, scenario_ids)
        Path(output).write_text(json.dumps(document))
        return document

    monkeypatch.setattr(ortholog_scope, "validate_scope_bundle", fake_validate)
    monkeypatch.setattr(
        ortholog_scope, "finalize_mapping_registry", fake_finalize,
    )
    monkeypatch.setattr(acquisition, "_safe_cache_root", lambda path: Path(path))
    monkeypatch.setattr(
        acquisition,
        "acquire",
        lambda *_args, **_kwargs: ({}, {"verified": True}),
    )

    def build_action(scenario_id):
        return [
            sys.executable,
            "-B",
            "-m",
            "benchmarks.compare.ortholog_scope",
            "build",
            "source.gff3",
            "source.fa",
            "truth.gff3",
            "target.fa",
            str(ortholog_root / scenario_id),
        ]

    scheduled = []
    run_count = 0

    def fake_run(actions, **_kwargs):
        nonlocal run_count
        run_count += 1
        selected = [Path(action[9]).name for action in actions]
        scheduled.append(selected)
        if run_count == 1:
            publish_bundle(selected[0])
            raise acquisition.AcquisitionError("simulated partial action failure")
        for scenario_id in selected:
            publish_bundle(scenario_id)
        return [
            {"index": index, "returncode": 0}
            for index, _scenario_id in enumerate(selected, start=1)
        ]

    monkeypatch.setattr(acquisition, "_run_bootstrap_actions", fake_run)

    def fake_materialize(*_args, ortholog_registry, **_kwargs):
        mapped = set()
        if ortholog_registry is not None:
            mapped = set(json.loads(
                Path(ortholog_registry).read_text()
            )["mappings"])
        missing = [
            scenario_id for scenario_id in scenario_ids
            if scenario_id not in mapped
        ]
        if missing:
            report = {
                "campaign_ready": False,
                "blockers": [
                    {"code": "missing_ortholog_map", "scenario_id": value}
                    for value in missing
                ],
                "remaining_actions": [
                    build_action(scenario_id) for scenario_id in missing
                ],
            }
        else:
            report = {
                "campaign_ready": True,
                "blockers": [],
                "remaining_actions": [],
                "registries": {
                    "benchmark": acquisition._file_record(benchmark_overlay),
                    "dataset": acquisition._file_record(dataset_overlay),
                },
            }
        preflight.parent.mkdir(parents=True, exist_ok=True)
        preflight.write_text(json.dumps(report))
        return report

    monkeypatch.setattr(acquisition, "materialize", fake_materialize)
    arguments = {
        "manifest_path": manifest_path,
        "cache_root": cache_root,
        "lock_path": cache_root / acquisition.LOCK_NAME,
        "benchmark_registry": tmp_path / "benchmarks.json",
        "dataset_registry": tmp_path / "datasets.json",
        "threads": 1,
        "max_active": 2,
        "accept_identity_pinned_bytes": True,
    }

    with pytest.raises(
        acquisition.AcquisitionError, match="simulated partial action failure",
    ):
        acquisition.bootstrap(manifest, **arguments)
    result = acquisition.bootstrap(manifest, **arguments)

    assert scheduled == [list(scenario_ids), ["bio-two"]]
    assert result["campaign_ready"] is True
    assert result["ortholog_bootstrap"]["valid_scenarios"] == list(
        scenario_ids
    )


def test_bootstrap_rejects_invalid_or_interrupted_bundle_with_quarantine(
        tmp_path, monkeypatch):
    from benchmarks.compare import ortholog_scope

    manifest = {
        "scenarios": [{"id": "bio-one", "kind": "biological"}],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    ortholog_root = tmp_path / "ortholog-scopes"
    invalid = ortholog_root / "bio-one"
    invalid.mkdir(parents=True)
    (invalid / "partial.txt").write_text("incomplete")
    monkeypatch.setattr(
        ortholog_scope,
        "validate_scope_bundle",
        lambda _path: (_ for _ in ()).throw(
            ortholog_scope.ScopeBuildError("missing manifest")
        ),
    )

    with pytest.raises(
        acquisition.AcquisitionError, match="quarantine.*bio-one",
    ):
        acquisition._bootstrap_ortholog_state(
            manifest,
            manifest_path=manifest_path,
            ortholog_root=ortholog_root,
            final_registry=ortholog_root / "ortholog_registry.json",
        )

    for path in invalid.iterdir():
        path.unlink()
    invalid.rmdir()
    interrupted = ortholog_root / ".bio-one.tmp-crash"
    interrupted.mkdir()
    with pytest.raises(
        acquisition.AcquisitionError, match="quarantine.*tmp-crash",
    ):
        acquisition._bootstrap_ortholog_state(
            manifest,
            manifest_path=manifest_path,
            ortholog_root=ortholog_root,
            final_registry=ortholog_root / "ortholog_registry.json",
        )


def test_bytecode_safe_action_command_is_idempotent():
    command = [sys.executable, "-B", "-m", "example"]
    assert acquisition._bytecode_safe_command(command) == command
    assert acquisition._bytecode_safe_command([
        sys.executable, "-m", "example",
    ])[1] == "-B"
