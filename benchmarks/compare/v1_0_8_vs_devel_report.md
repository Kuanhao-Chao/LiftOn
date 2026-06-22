# LiftOn v1.0.8 (previous stable release) vs `devel` HEAD — comprehensive comparison

**stable** = the real LiftOn **v1.0.8** binary (git tag `v1.0.8` = `e503643d`, what `pip install lifton` delivers), built into an isolated `lifton_stable` conda env. **devel** = current `devel` HEAD (Iteration 19, `4496dd5`), **41 commits ahead**. **devel_legacy** = the devel binary run with `--legacy-merge --full-dp-align --gene-only` (emulates v1.0.8's *algorithm* defaults) — used only on the controlled arm to attribute deltas.

Accuracy is the **version-agnostic neutral re-score** (`benchmarks.compare.evaluator`): each predicted protein is re-aligned to the reference protein with LiftOn's own parasail kernel, identically for both versions — it never trusts either tool's self-reported `protein_identity`. Validity uses the **devel** `gff3-validate` on both outputs (one consistent yardstick). `-V` is NOT a version discriminator (devel's `__version__` still reads `v1.0.8`); provenance is pinned on each install's git SHA + `lifton.__file__`.

## Executive summary

Controlled arm (cached `-L`/`-M`, `-t 1` both — isolates LiftOn's in-process logic):

| dataset | neutral mean PI stable→devel (Δ) | impr/regr | recovered +dev/−dev | wall stable→devel | RSS stable→devel |
|---|---|---|---|---|---|
| drosophila | 0.9221 → 0.92578 (+0.003680) | 115/1 | +0/−0 | 291.9→220.2s | 2364.5→536.4 MB |
| mouse_to_rat | 0.90521 → 0.9092 (+0.003990) | 83/0 | +0/−0 | 268→131.9s | 10690.6→538.3 MB |
| rice | 0.99905 → 0.99905 (+0.000000) | 0/0 | +0/−0 | 151.2→135.8s | 1399.9→508.3 MB |
| arabidopsis | 0.99886 → 0.99897 (+0.000110) | 9/0 | +0/−0 | 266.3→256.9s | 1367.8→629.5 MB |
| bee | 0.9916 → 0.9916 (+0.000000) | 0/0 | +0/−0 | 230.5→107.1s | 10939.0→457.3 MB |

## Arm: controlled

### Accuracy (neutral re-score)

| dataset | n coding | mean PI stable | mean PI devel | Δ mean | median s→d | %identical s→d | improved | regressed | net/txn | self−neutral bias (s/d) |
|---|---|---|---|---|---|---|---|---|---|---|
| drosophila | 7251 | 0.9221 | 0.92578 | +0.003680 | 0.95793→0.958 | 0.04758→0.04758 | 115 | 1 | +0.003677 | +0.0001/+0.0000 |
| mouse_to_rat | 2343 | 0.90521 | 0.9092 | +0.003990 | 0.95946→0.9596 | 0.05705→0.05705 | 83 | 0 | +0.003996 | +0.0000/+0.0000 |
| rice | 5850 | 0.99905 | 0.99905 | +0.000000 | 1→1 | 0.99538→0.99538 | 0 | 0 | +0.000000 | -0.0000/-0.0000 |
| arabidopsis | 12653 | 0.99886 | 0.99897 | +0.000110 | 1→1 | 0.99533→0.99533 | 9 | 0 | +0.000115 | -0.0000/-0.0000 |
| bee | 3130 | 0.9916 | 0.9916 | +0.000000 | 1→1 | 0.73062→0.73062 | 0 | 0 | +0.000000 | +0.0000/+0.0000 |

### Completeness

| dataset | recovered coding stable→devel | completeness_coding s→d | feature_total s→d | recovered +dev/−dev |
|---|---|---|---|---|
| drosophila | 7131→7131 | 0.98345→0.98345 | 0.91009→0.91009 | +0/−0 |
| mouse_to_rat | 2126→2126 | 0.90738→0.90738 | 0.74875→0.74875 | +0/−0 |
| rice | 5850→5850 | 1→1 | 0.97572→0.97572 | +0/−0 |
| arabidopsis | 12632→12632 | 0.99834→0.99834 | 0.91601→0.91601 | +0/−0 |
| bee | 3122→3122 | 0.99744→0.99744 | 0.96665→0.96665 | +0/−0 |

### Performance

| dataset | wall stable | wall devel | wall speedup | RSS stable | RSS devel | CPU stable | CPU devel |
|---|---|---|---|---|---|---|---|
| drosophila | 291.85s | 220.25s | 1.33x | 2364.5 MB | 536.4 MB | 275.2s | 203.9s |
| mouse_to_rat | 268.01s | 131.94s | 2.03x | 10690.6 MB | 538.3 MB | 264.6s | 128.5s |
| rice | 151.22s | 135.8s | 1.11x | 1399.9 MB | 508.3 MB | 139.3s | 124.0s |
| arabidopsis | 266.34s | 256.9s | 1.04x | 1367.8 MB | 629.5 MB | 238.3s | 224.1s |
| bee | 230.48s | 107.1s | 2.15x | 10939.0 MB | 457.3 MB | 227.7s | 104.8s |

### Validity (devel `gff3-validate`, one yardstick on both outputs)

| dataset | stable errors/warnings | devel errors/warnings |
|---|---|---|
| drosophila | 52/78 | 51/63 |
| mouse_to_rat | 50/8 | 50/5 |
| rice | 50/10 | 50/10 |
| arabidopsis | 50/18 | 50/18 |
| bee | 50/0 | 50/0 |

### Attribution (controlled): devel_legacy emulates v1.0.8's algorithm defaults

`devel_legacy` = devel + `--legacy-merge --full-dp-align --gene-only`. **devel_legacy vs stable** = the residual from the NON-flag-gated changes (Phase-4 gene-ID fix + Phase-13.5C canonical writer). **devel vs devel_legacy** = the flag-gated algorithm promotions (best-of-outcome merge + band-everything alignment + gene-like lift).

| dataset | mean PI stable | mean PI devel_legacy | mean PI devel | residual (legacy−stable) impr/regr | promotions (devel−legacy) impr/regr |
|---|---|---|---|---|---|
| drosophila | 0.9221 | 0.92208 | 0.92578 | 0/2 (net -0.000063) | 115/1 (net +0.003741) |
| mouse_to_rat | 0.90521 | 0.90398 | 0.9092 | 0/12 (net -0.001225) | 91/0 (net +0.005221) |
| rice | 0.99905 | 0.99905 | 0.99905 | 0/0 (net +0.000000) | 0/0 (net +0.000000) |
| arabidopsis | 0.99886 | 0.99886 | 0.99897 | 0/0 (net +0.000000) | 9/0 (net +0.000115) |
| bee | 0.9916 | 0.9916 | 0.9916 | 0/0 (net +0.000000) | 0/0 (net +0.000000) |

## Arm: fresh

### Accuracy (neutral re-score)

| dataset | n coding | mean PI stable | mean PI devel | Δ mean | median s→d | %identical s→d | improved | regressed | net/txn | self−neutral bias (s/d) |
|---|---|---|---|---|---|---|---|---|---|---|
| drosophila | 7251 | 0.92567 | 0.92693 | +0.001260 | 0.9585→0.95873 | 0.04898→0.04899 | 60 | 0 | +0.001260 | +0.0001/+0.0000 |
| arabidopsis | 12653 | 0.99893 | 0.999 | +0.000070 | 1→1 | 0.99541→0.99541 | 3 | 0 | +0.000071 | -0.0000/-0.0000 |

### Completeness

| dataset | recovered coding stable→devel | completeness_coding s→d | feature_total s→d | recovered +dev/−dev |
|---|---|---|---|---|
| drosophila | 7131→7130 | 0.98345→0.98331 | 0.91009→0.91161 | +1/−2 |
| arabidopsis | 12632→12633 | 0.99834→0.99842 | 0.91601→0.9888 | +659/−1 |

Per-type reference-feature recovery (stable→devel `n_recovered`), types where they differ:

| dataset | feature type | n_reference | stable recovered | devel recovered |
|---|---|---|---|---|
| drosophila | CDS | 37963 | 1761 | 2050 |
| drosophila | exon | 44245 | 19932 | 20168 |
| drosophila | mRNA | 7251 | 7131 | 7130 |
| drosophila | pseudogene | 44 | 0 | 21 |
| arabidopsis | CDS | 77644 | 12575 | 12576 |
| arabidopsis | exon | 87301 | 85056 | 86768 |
| arabidopsis | mRNA | 13336 | 12632 | 13290 |
| arabidopsis | pseudogene | 930 | 0 | 903 |
| arabidopsis | transcript | 450 | 170 | 449 |

### Performance

| dataset | wall stable | wall devel | wall speedup | RSS stable | RSS devel | CPU stable | CPU devel |
|---|---|---|---|---|---|---|---|
| drosophila | 414.33s | 320.58s | 1.29x | 2558.7 MB | 1208.9 MB | 428.0s | 346.0s |
| arabidopsis | 336.82s | 297.02s | 1.13x | 2325.2 MB | 2436.0 MB | 346.9s | 321.4s |

### Validity (devel `gff3-validate`, one yardstick on both outputs)

| dataset | stable errors/warnings | devel errors/warnings |
|---|---|---|
| drosophila | 126/63 | 75/64 |
| arabidopsis | 100/18 | 50/77 |

## Arm: full

> **Recorded crashes:** `rice` — stable crashed; `arabidopsis` — stable crashed. A crashed version's cells show `n/a`; the surviving version still scored. (v1.0.8 crashes on some full genomes with an unhandled `gffutils.FeatureNotFoundError` whose error handler itself raises `TypeError: __str__ returned non-string` — a robustness bug devel fixed.)

### Accuracy (neutral re-score)

| dataset | n coding | mean PI stable | mean PI devel | Δ mean | median s→d | %identical s→d | improved | regressed | net/txn | self−neutral bias (s/d) |
|---|---|---|---|---|---|---|---|---|---|---|
| rice | 42580 | n/a | 0.99818 | n/a | n/a→1 | n/a→0.98791 | ? | ? | n/a | n/a/-0.0000 |
| arabidopsis | 48265 | n/a | 0.99904 | n/a | n/a→1 | n/a→0.99369 | ? | ? | n/a | n/a/-0.0000 |
| bee | 23471 | 0.99024 | 0.99033 | +0.000090 | 1→1 | 0.70283→0.70275 | 55 | 8 | +0.000088 | +0.0000/-0.0000 |

### Completeness

| dataset | recovered coding stable→devel | completeness_coding s→d | feature_total s→d | recovered +dev/−dev |
|---|---|---|---|---|
| rice | None→42527 | n/a→0.99876 | n/a→0.97587 | +?/−? |
| arabidopsis | None→48207 | n/a→0.9988 | n/a→0.97951 | +?/−? |
| bee | 23364→23364 | 0.99544→0.99544 | 0.94493→0.94565 | +0/−0 |

Per-type reference-feature recovery (stable→devel `n_recovered`), types where they differ:

| dataset | feature type | n_reference | stable recovered | devel recovered |
|---|---|---|---|---|
| rice | CDS | 238915 | 0 | 42032 |
| rice | exon | 334017 | 0 | 331837 |
| rice | gene | 33624 | 0 | 33565 |
| rice | lnc_RNA | 6277 | 0 | 6271 |
| rice | mRNA | 42580 | 0 | 42527 |
| rice | pseudogene | 1605 | 0 | 1593 |
| rice | pseudogenic_tRNA | 3 | 0 | 2 |
| rice | rRNA | 200 | 0 | 188 |
| rice | snRNA | 71 | 0 | 71 |
| rice | snoRNA | 610 | 0 | 610 |
| rice | tRNA | 672 | 0 | 666 |
| rice | transcript | 3840 | 0 | 3836 |
| arabidopsis | CDS | 286264 | 0 | 284355 |
| arabidopsis | antisense_RNA | 92 | 0 | 92 |
| arabidopsis | exon | 324848 | 0 | 322744 |
| arabidopsis | gene | 33468 | 0 | 33175 |
| arabidopsis | lnc_RNA | 3878 | 0 | 3878 |
| arabidopsis | mRNA | 52177 | 0 | 52086 |
| arabidopsis | miRNA | 428 | 0 | 1 |
| arabidopsis | ncRNA | 540 | 0 | 314 |
| arabidopsis | primary_transcript | 326 | 0 | 325 |
| arabidopsis | pseudogene | 4851 | 0 | 4816 |
| arabidopsis | rRNA | 14 | 0 | 14 |
| arabidopsis | snRNA | 82 | 0 | 77 |
| arabidopsis | snoRNA | 287 | 0 | 283 |
| arabidopsis | tRNA | 684 | 0 | 663 |
| arabidopsis | three_prime_UTR | 23 | 0 | 1 |
| arabidopsis | transcript | 1826 | 0 | 1822 |
| bee | CDS | 226465 | 16432 | 138297 |
| bee | exon | 267977 | 189742 | 194744 |
| bee | gene | 12356 | 12296 | 12295 |
| bee | pseudogene | 42 | 0 | 33 |
| bee | rRNA | 57 | 56 | 55 |

### Performance

| dataset | wall stable | wall devel | wall speedup | RSS stable | RSS devel | CPU stable | CPU devel |
|---|---|---|---|---|---|---|---|
| rice | Nones | 9807.0s | n/a | None MB | 8150.3 MB | Nones | 26173.9s |
| arabidopsis | Nones | 11772.0s | n/a | None MB | 12556.1 MB | Nones | 30727.6s |
| bee | 3395.18s | 7504.0s | 0.45x | 23083.1 MB | 7098.5 MB | 3558.0s | 17238.2s |

### Validity (devel `gff3-validate`, one yardstick on both outputs)

| dataset | stable errors/warnings | devel errors/warnings |
|---|---|---|
| rice | ?/? | 129/102 |
| arabidopsis | ?/? | 107/114 |
| bee | 103/96 | 53/99 |

