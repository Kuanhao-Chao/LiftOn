"""Validate canonical-v2 dataset provenance and acquisition locking."""
from __future__ import annotations

import copy
import hashlib
from pathlib import Path

import pytest

from benchmarks import manifest_tools


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = (
    REPO_ROOT / "benchmarks" / "manifests" / "canonical_v2_datasets.json"
)


@pytest.fixture(scope="module")
def manifest():
    return manifest_tools.validate_manifest(
        manifest_tools.load_json(MANIFEST_PATH)
    )


def test_manifest_pins_the_twelve_approved_scenarios(manifest):
    assert tuple(row["id"] for row in manifest["scenarios"]) == (
        manifest_tools.EXPECTED_SCENARIO_IDS
    )
    assert len(manifest["sources"]) == 18
    assert [source["id"] for source in manifest["sources"]] == sorted(
        source["id"] for source in manifest["sources"]
    )
    assert {row["kind"] for row in manifest["scenarios"]} == {
        "biological", "synthetic", "protocol",
    }


def test_manifest_uses_fixed_source_identities_and_no_fake_hashes(manifest):
    remote = [
        source for source in manifest["sources"]
        if source["transport"] != "existing_registry"
    ]
    assert len(remote) == 15
    for source in remote:
        for file_record in source["files"]:
            assert file_record["expected_sha256"] is None
            assert file_record["pin_state"] == "identity_pinned_bytes_pending"
            if source["transport"] == "https_files":
                assert file_record["locator"]["url"].startswith("https://")
                assert "?" not in file_record["locator"]["url"]
    by_id = {source["id"]: source for source in manifest["sources"]}
    assert by_id["ncbi_grch38_p14"]["identity"]["accession"] == "GCF_000001405.40"
    assert by_id["ncbi_chm13_v2"]["identity"]["accession"] == "GCF_009914755.1"
    assert by_id["ensembl_grch38_release116"]["identity"]["release"] == 116
    assert {
        row["role"]: row["locator"]["url"]
        for row in by_id["ensembl_grch38_release116"]["files"]
    } == {
        "annotation_gtf": (
            "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh38.116.gtf.gz"
        ),
        "genome": (
            "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
    }
    assert by_id["flybase_dmel_fb2026_02"]["release"] == "FB2026_02 r6.68"
    assert by_id["wormbase_celegans_wbps19"]["release"].startswith("WBPS19")

    scenarios = {row["id"]: row for row in manifest["scenarios"]}
    assert scenarios["v2_dialect_ensembl116_gtf"]["inputs"][
        "reference_genome"
    ] == "ensembl_grch38_release116:genome"


def test_synthetic_designs_are_exact_builder_inputs(manifest):
    scenarios = {row["id"]: row for row in manifest["scenarios"]}
    fragmented = scenarios["v2_synth_chr22_fragmented"]
    sv = scenarios["v2_synth_chr22_sv"]

    for scenario in (fragmented, sv):
        assert set(scenario["inputs"]) == {
            "reference_annotation", "reference_genome",
        }
        assert "base_target_genome" not in scenario["inputs"]

    assert "seed" not in fragmented["design"]
    assert manifest_tools.synthetic_builder_kwargs(fragmented) == {
        "source_seqid": "chr22",
        "cuts": (
            17_290_898, 22_364_968, 27_189_473, 34_494_491,
            41_032_950, 45_086_965, 47_832_626, 49_262_605,
        ),
    }
    assert manifest_tools.expected_synthetic_transform(fragmented) == {
        "id": "v2_synth_chr22_fragmented",
        "kind": "fragmentation",
        "cuts_after_source_coordinate": [
            17_290_898, 22_364_968, 27_189_473, 34_494_491,
            41_032_950, 45_086_965, 47_832_626, 49_262_605,
        ],
    }

    assert manifest_tools.synthetic_builder_kwargs(sv) == {
        "source_seqid": "chr22",
        "deletion": (18_177_398, 18_516_337),
        "inversion": (36_140_323, 38_318_084),
        "duplication": (41_092_592, 44_331_714),
        "insert_after": 45_086_965,
        "insertion_length": 100_000,
        "seed": "lifton-canonical-v2-20260718",
    }
    assert manifest_tools.expected_synthetic_transform(sv)[
        "insertion_sha256"
    ] == "aa96e5dc55452f08fec81cee74cd299c42ad8187d3c472bd240a15d470eb5236"


def test_manifest_rejects_noncanonical_synthetic_operations(manifest):
    arbitrary = copy.deepcopy(manifest)
    arbitrary["scenarios"][8]["design"]["operations"] = [{"type": "bogus"}]
    with pytest.raises(manifest_tools.ManifestError, match="keys must be exactly"):
        manifest_tools.validate_manifest(arbitrary)

    invalid_cuts = copy.deepcopy(manifest)
    invalid_cuts["scenarios"][8]["design"]["operations"][0]["cut_after"][1] = (
        invalid_cuts["scenarios"][8]["design"]["operations"][0]["cut_after"][0]
    )
    with pytest.raises(manifest_tools.ManifestError, match="strictly increasing"):
        manifest_tools.validate_manifest(invalid_cuts)

    missing_sv_parameter = copy.deepcopy(manifest)
    del missing_sv_parameter["scenarios"][9]["design"]["operations"][3][
        "generator"
    ]
    with pytest.raises(manifest_tools.ManifestError, match="keys must be exactly"):
        manifest_tools.validate_manifest(missing_sv_parameter)

    changed_seed = copy.deepcopy(manifest)
    changed_seed["scenarios"][9]["design"]["seed"] = "unsafe seed"
    with pytest.raises(manifest_tools.ManifestError, match="safe non-empty"):
        manifest_tools.validate_manifest(changed_seed)


def test_synthetic_builder_follows_valid_operation_and_seed_changes(manifest):
    changed = copy.deepcopy(manifest)
    changed["scenarios"][8]["design"]["operations"][0]["cut_after"][0] += 1
    changed["scenarios"][9]["design"]["operations"][0]["start"] += 1
    changed["scenarios"][9]["design"]["seed"] += "-alternate"
    validated = manifest_tools.validate_manifest(changed)
    fragmented = validated["scenarios"][8]
    sv = validated["scenarios"][9]

    assert manifest_tools.synthetic_builder_kwargs(fragmented)["cuts"][0] == (
        17_290_899
    )
    assert manifest_tools.expected_synthetic_transform(fragmented)[
        "cuts_after_source_coordinate"
    ][0] == 17_290_899

    kwargs = manifest_tools.synthetic_builder_kwargs(sv)
    assert kwargs["deletion"] == (18_177_399, 18_516_337)
    assert kwargs["seed"] == "lifton-canonical-v2-20260718-alternate"
    transform = manifest_tools.expected_synthetic_transform(sv)
    assert transform["deletion"] == [18_177_399, 18_516_337]
    assert transform["deleted_length"] == 338_939
    assert transform["seed"] == "lifton-canonical-v2-20260718-alternate"
    assert transform["insertion_sha256"] != (
        "aa96e5dc55452f08fec81cee74cd299c42ad8187d3c472bd240a15d470eb5236"
    )


def test_manifest_rejects_synthetic_target_binding(manifest):
    changed = copy.deepcopy(manifest)
    changed["scenarios"][8]["inputs"]["base_target_genome"] = (
        "existing_human_mane:target_genome"
    )
    with pytest.raises(
        manifest_tools.ManifestError,
        match="synthetic inputs must be exactly",
    ):
        manifest_tools.validate_manifest(changed)


def test_acquisition_dry_run_is_deterministic_and_non_mutating(manifest, tmp_path):
    cache_root = tmp_path / "must-not-be-created"
    first = manifest_tools.build_acquisition_plan(manifest, cache_root)
    second = manifest_tools.build_acquisition_plan(manifest, cache_root)
    assert first == second
    assert first["dry_run"] is True
    assert first["remote_source_count"] == 15
    assert len(first["steps"]) == 15
    assert not cache_root.exists()
    assert all(step["commands"] for step in first["steps"])
    assert {
        step["transport"] for step in first["steps"]
    } == {"ncbi_datasets", "https_files"}
    assert len(first["plan_sha256"]) == 64


def test_manifest_rejects_changed_scenario_order_and_unknown_binding(manifest):
    reordered = copy.deepcopy(manifest)
    reordered["scenarios"][0], reordered["scenarios"][1] = (
        reordered["scenarios"][1], reordered["scenarios"][0]
    )
    with pytest.raises(manifest_tools.ManifestError, match="IDs/order"):
        manifest_tools.validate_manifest(reordered)

    unknown = copy.deepcopy(manifest)
    unknown["scenarios"][0]["inputs"]["target_truth"] = "missing:annotation"
    with pytest.raises(manifest_tools.ManifestError, match="unknown binding"):
        manifest_tools.validate_manifest(unknown)


def _create_lock(manifest, cache_root):
    sources = {}
    for source in manifest["sources"]:
        if source["transport"] == "existing_registry":
            continue
        files = {}
        for file_record in source["files"]:
            role = file_record["role"]
            payload = f"{source['id']}:{role}\n".encode()
            digest = hashlib.sha256(payload).hexdigest()
            relative = Path("sha256") / digest[:2] / digest / f"{source['id']}-{role}"
            path = cache_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            files[role] = {
                "path": str(relative),
                "bytes": len(payload),
                "sha256": digest,
            }
        sources[source["id"]] = {
            "request_sha256": manifest_tools.source_request_sha256(source),
            "files": files,
        }
    return {
        "schema_version": 1,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "sources": sources,
    }


def test_lock_verifier_rehashes_content_addressed_files(manifest, tmp_path):
    lock = _create_lock(manifest, tmp_path)
    result = manifest_tools.verify_acquisition_lock(manifest, lock, tmp_path)
    assert result["verified"] is True
    assert result["source_count"] == 15
    assert result["file_count"] == sum(
        len(source["files"]) for source in manifest["sources"]
        if source["transport"] != "existing_registry"
    )

    source_id = sorted(lock["sources"])[0]
    role = sorted(lock["sources"][source_id]["files"])[0]
    path = tmp_path / lock["sources"][source_id]["files"][role]["path"]
    path.write_bytes(b"tampered\n")
    with pytest.raises(
        manifest_tools.ManifestError,
        match="byte count changed|content changed",
    ):
        manifest_tools.verify_acquisition_lock(manifest, lock, tmp_path)


def test_manifest_rejects_unpinned_identity_without_pending_policy(manifest):
    changed = copy.deepcopy(manifest)
    changed["sources"][0]["files"][0]["pin_state"] = "sha256_pinned"
    with pytest.raises(manifest_tools.ManifestError, match="require lowercase SHA"):
        manifest_tools.validate_manifest(changed)
