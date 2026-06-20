"""Validate the tiered benchmark registry (``benchmarks/tiers.json``) and the
``register_tiers`` merge logic. Pure-Python — no network, no NCBI ``datasets``
CLI, no heavy deps — so it runs in the normal pytest env alongside the rest of
the suite.
"""
import json
from pathlib import Path

import pytest

from benchmarks.register_tiers import build_entries, build_entry, merge

REPO = Path(__file__).resolve().parents[1]
TIERS_JSON = REPO / "benchmarks" / "tiers.json"
BENCHMARKS_JSON = REPO / "benchmarks" / "compare" / "benchmarks.json"

# The four tiers reuse the repo's existing divergence_class vocabulary.
VALID_TIERS = {
    "same_species",
    "close_cross_species",
    "distant_cross_species",
    "very_distant_cross_species",
}
STD_ENTRY_KEYS = {
    "id", "species", "cross_species", "dimension", "divergence_class",
    "ref_genome", "ref_gff", "ref_proteins", "tgt_genome",
    "ref_chrom", "tgt_chrom", "miniprot_target_space", "annotation_database",
}


@pytest.fixture(scope="module")
def reg():
    with open(TIERS_JSON) as fh:
        return json.load(fh)


def test_registry_parses_and_has_12_pairs(reg):
    assert isinstance(reg["pairs"], list)
    assert len(reg["pairs"]) == 12


def test_tiers_map_matches_vocabulary(reg):
    assert set(reg["tiers"]) == VALID_TIERS
    # 1..4 numbering, one per tier
    assert sorted(t["num"] for t in reg["tiers"].values()) == [1, 2, 3, 4]


def test_three_pairs_per_tier(reg):
    from collections import Counter
    counts = Counter(p["tier"] for p in reg["pairs"])
    assert set(counts) == VALID_TIERS
    assert all(n == 3 for n in counts.values()), counts


def test_pair_ids_unique(reg):
    ids = [p["id"] for p in reg["pairs"]]
    assert len(ids) == len(set(ids))


@pytest.mark.parametrize("field", [
    "id", "tier", "dimension", "ref_species", "tgt_species",
    "cross_species", "divergence_mya", "justification",
    "ref", "tgt", "ref_chrom", "tgt_chrom",
])
def test_every_pair_has_required_fields(reg, field):
    for p in reg["pairs"]:
        assert field in p, f"{p.get('id')} missing {field}"


def test_each_pair_has_accession_or_local_hint(reg):
    for p in reg["pairs"]:
        ref, tgt = p["ref"], p["tgt"]
        # an accession is always present (the canonical / download fallback)
        assert ref.get("acc"), f"{p['id']} ref has no accession"
        assert tgt.get("acc"), f"{p['id']} tgt has no accession"
        # local hints, when present, are absolute paths
        for hint in ("local_fna", "local_gff", "local_faa"):
            v = ref.get(hint)
            if v:
                assert str(v).startswith("/"), f"{p['id']} ref.{hint} not absolute"
        if tgt.get("local_fna"):
            assert str(tgt["local_fna"]).startswith("/")


def test_tier_consistency_with_divergence(reg):
    # the new pairs fall cleanly inside their MYA bands
    for p in reg["pairs"]:
        mya = p["divergence_mya"]
        tier = p["tier"]
        if tier == "same_species":
            assert mya == 0 and p["cross_species"] is False
        elif tier == "close_cross_species":
            assert 0 < mya < 20 and p["cross_species"] is True
        elif tier == "distant_cross_species":
            assert 20 <= mya <= 60 and p["cross_species"] is True
        elif tier == "very_distant_cross_species":
            assert mya > 60 and p["cross_species"] is True


def test_very_distant_targets_use_whole(reg):
    # synteny is gone past ~60 MYA -> the syntenic-target picker is skipped
    for p in reg["pairs"]:
        if p["tier"] == "very_distant_cross_species":
            assert p["tgt_chrom"] == "WHOLE", p["id"]


def test_build_entry_has_standard_schema(reg):
    for p in reg["pairs"]:
        e = build_entry(p, reg["data_root"])
        assert set(e) == STD_ENTRY_KEYS, p["id"]
        assert e["divergence_class"] in VALID_TIERS
        assert e["ref_genome"].endswith("/ref.fna")
        assert e["ref_gff"].endswith("/ref.gff")
        assert e["tgt_genome"].endswith("/tgt.fna")
        assert reg["data_root"] in e["ref_genome"]
        assert p["tier"] in e["ref_genome"]


def test_build_entries_count(reg):
    assert len(build_entries(reg)) == 12


def test_merge_is_idempotent_and_collision_free(reg):
    entries = build_entries(reg)
    existing = [{"id": "human_mane"}, {"id": "drosophila"}]  # stand-in registry
    once = merge(existing, entries)
    twice = merge(once, entries)
    assert once == twice                       # idempotent
    ids = [e["id"] for e in twice]
    assert len(ids) == len(set(ids))           # no duplicates
    # existing kept, in order, ahead of the appended tier entries
    assert ids[:2] == ["human_mane", "drosophila"]
    assert ids[2:] == [e["id"] for e in entries]


def test_existing_benchmarks_json_still_parses():
    # guard: the committed registry the merge targets is valid JSON
    with open(BENCHMARKS_JSON) as fh:
        bench = json.load(fh)
    assert isinstance(bench["benchmarks"], list)
