# LiftOn benchmark expansion — 4-tier species-pair set

A curated expansion of the LiftOn full-genome annotation-liftover benchmark with
**12 reference→target genome pairs across four evolutionary-distance tiers**.
Each pair is a *practical* liftover use case (a highly-annotated model/reference
genome → a related, less-annotated genome) and reuses the repo's existing
`divergence_class` vocabulary so it drops straight into the 4-way comparison
harness (`benchmarks/compare/benchmarks.json` → Liftoff / miniprot /
LiftOn-stable / LiftOn-devel, neutral parasail re-scorer).

The machine-readable source of truth is **`benchmarks/tiers.json`**; this file
is the human-readable companion. Accessions were verified against NCBI Datasets
and divergence times against TimeTree (2026-06-18).

## 1. Distance tiers

| Tier | Name | MYA band | `divergence_class` | What it stresses |
|---|---|---|---|---|
| 1 | Intra-species / strain | ~0 | `same_species` | Cultivar/strain-level transfer; near-perfect synteny, SV-driven differences. |
| 2 | Close cross-species | < 20 | `close_cross_species` | Recently diverged relatives; high identity, occasional gene gain/loss. |
| 3 | Moderate cross-species | 20–60 | `distant_cross_species` | Mid-range divergence; where protein-aware chaining starts to pay off. |
| 4 | Distant cross-species | > 60 | `very_distant_cross_species` | Deep divergence; DNA synteny largely gone (`tgt_chrom=WHOLE`), protein signal dominates. |

The tier directory names equal the `divergence_class` strings. The MYA bands map
the user's four tiers onto the registry's existing four classes; existing pairs
are **not** relabelled.

## 2. Divergence summary table

`local` in the source column = the genome is already on disk and is symlinked
(not re-downloaded); the listed accession is the canonical identity / download
fallback. Reference accessions carry a RefSeq GFF3; targets need only the genome.

### Tier 1 — intra-species / strain (`same_species`, ~0 MYA)
| Pair `id` | Reference (acc · assembly) | Target (acc · assembly) | MYA | Why it's useful |
|---|---|---|---|---|
| `t1_maize_b73_to_mo17` | *Zea mays* B73 · `GCF_902167145.1` (Zm-B73-NAM-5.0) | *Zea mays* Mo17 · `GCA_022117705.1` (Zm-Mo17 T2T) | ~0 | Lift the gold-standard B73 NAM annotation onto the other foundational maize inbred; SV-rich crop pair in daily NAM/pan-genome use. |
| `t1_soybean_w82_to_lee` | *Glycine max* Williams 82 · `GCF_000004515.6` (v4.0) | *Glycine max* Lee · `GCA_002905335.2` (glyma.Lee.gnm1 v2) | ~0 | Reference→cultivar transfer in the major legume crop; soybean pan-genome work. |
| `t1_tomato_microtom_to_heinz` | *Solanum lycopersicum* Micro-Tom · `GCF_036512215.1` (SLM_r2.1) | *S. lycopersicum* Heinz 1706 · `GCA_000188115.5` (SL4.0) | ~0 | Cross-cultivar transfer in the model Solanaceae (current RefSeq is Micro-Tom; Heinz is the classic genome). |

### Tier 2 — close cross-species (`close_cross_species`, < 20 MYA)
| Pair `id` | Reference (acc · assembly) | Target (acc · assembly) | MYA | Why it's useful |
|---|---|---|---|---|
| `t2_human_to_gorilla` | *Homo sapiens* GRCh38 · `GCF_000001405.40` (local) | *Gorilla gorilla* · `GCF_029281585.2` (NHGRI_mGorGor1-v2.1, T2T) | ~8.6 | Human RefSeq → T2T great-ape; primate comparative / conservation genomics. |
| `t2_mouse_to_caroli` | *Mus musculus* GRCm39 · `GCF_000001635.27` (local) | *Mus caroli* · `GCF_900094665.2` (CAROLI_EIJ_v1.1) | ~4 | Mouse reference → close wild outgroup used in allele-specific / imprinting F1 studies. |
| `t2_tomato_to_potato` | *Solanum lycopersicum* Micro-Tom · `GCF_036512215.1` | *Solanum tuberosum* · `GCF_000226075.1` (SolTub_3.0) | ~8 | The canonical Solanaceae close-relative liftover for crop genomics. |

### Tier 3 — moderate cross-species (`distant_cross_species`, 20–60 MYA)
| Pair `id` | Reference (acc · assembly) | Target (acc · assembly) | MYA | Why it's useful |
|---|---|---|---|---|
| `t3_human_to_macaque` | *Homo sapiens* GRCh38 · `GCF_000001405.40` (local) | *Macaca mulatta* · `GCF_003339765.1` (Mmul_10) | ~29 | Human annotation → the dominant non-human-primate biomedical model. |
| `t3_dog_to_cat` | *Canis lupus familiaris* · `GCF_000002285.5` (Dog10K_Boxer_Tasha) | *Felis catus* · `GCF_018350175.1` (F.catus_Fca126_mat1.0) | ~50 | Carnivora / veterinary comparative genomics. |
| `t3_human_to_marmoset` | *Homo sapiens* GRCh38 · `GCF_000001405.40` (local) | *Callithrix jacchus* · `GCF_049354715.1` (calJac240_pri) | ~43 | Human → the major New-World-monkey neuroscience model (reference superseded late-2025). |

### Tier 4 — distant cross-species (`very_distant_cross_species`, > 60 MYA)
| Pair `id` | Reference (acc · assembly) | Target (acc · assembly) | MYA | Why it's useful |
|---|---|---|---|---|
| `t4_human_to_chicken` | *Homo sapiens* GRCh38 · `GCF_000001405.40` (local) | *Gallus gallus* · `GCF_016699485.2` (GRCg7b, local) | ~320 | Classic amniote deep-conservation stress test; both genomes already on disk. |
| `t4_human_to_xenopus` | *Homo sapiens* GRCh38 · `GCF_000001405.40` (local) | *Xenopus tropicalis* · `GCF_000004195.4` (UCB_Xtro_10.0) | ~350 | Human → frog developmental-biology model; tetrapod-scale divergence. |
| `t4_drosophila_to_bee` | *Drosophila melanogaster* · `GCF_000001215.4` (local) | *Apis mellifera* · `GCF_003254395.2` (Amel_HAv3.1, local) | ~345 | Holometabolous-insect comparative across Diptera↔Hymenoptera. |

The five **human** pairs pin `ref_chrom=chr20`, extending the registry's existing
controlled human ladder (chr20 → mouse / zebrafish), so the *same* human chr20
gene set is lifted across gorilla → macaque → marmoset → chicken → Xenopus
(plus the existing mouse/zebrafish) at increasing divergence. Other pairs use
`ref_chrom=AUTO_LARGEST_CODING`. Tiers 1–3 use `tgt_chrom=AUTO_SYNTENIC`
(minimap2 finds the syntenic target chromosome); tier 4 uses `tgt_chrom=WHOLE`
because genome-wide DNA synteny is gone at that distance.

*Verified-but-not-included extras* (easy to add to `tiers.json` if more breadth
is wanted): chicken→zebra-finch `GCF_003957565.2` (~100 MYA, avian) for tier 4;
pig `GCF_000003025.6`→cattle `GCF_002263795.3` (~62 MYA) at the tier-3/4 boundary.

## 3. Fetch & register

All paths land under `/ccb/salz3/kh.chao/lifton_benchmark_data/tiers/<tier>/<id>/`
with canonical names `ref.fna` / `ref.gff` / `ref.faa` (optional) / `tgt.fna`.
Run from the repo root with the pinned env:

```bash
export PYTHONNOUSERSITE=1
PY=/home/kh.chao/miniconda3/envs/lifton_devel/bin/python

# 1. download / symlink the inputs (idempotent, resumable; writes a manifest.json)
$PY -m benchmarks.fetch_tiered                      # all 12 pairs
$PY -m benchmarks.fetch_tiered t2_human_to_gorilla  # one pair
$PY -m benchmarks.fetch_tiered close_cross_species  # one whole tier

# 2. register the pairs into the 4-way comparison registry
$PY -m benchmarks.register_tiers --dry-run          # preview
$PY -m benchmarks.register_tiers                     # write benchmarks.json
```

`fetch_tiered.py` reuses the proven `benchmarks/compare/fetch_new_pairs.py`
download helpers (NCBI `datasets` CLI, 3× retry, idempotent extract). References
already on disk (human, mouse, drosophila) and the chicken/bee genomes are
symlinked rather than re-downloaded; each unique accession is cached so a shared
genome (the tomato Micro-Tom reference) downloads once. Per-pair failures are
logged and skipped — they do not abort the run — and `manifest.json` records the
source / byte-size / status of every canonical file.

Once registered, each pair is run through the existing 4-way comparison harness
(`benchmarks/compare/run_compare.py` / `fourway_compare.py`) like any other cell.
