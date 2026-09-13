# Release readiness — v1.0.12

Written 2026-08-25. This release responds to [GitHub issue
#71](https://github.com/Kuanhao-Chao/LiftOn/issues/71) with a resource-aware
native-aligner schedule, bounded failure diagnostics, and a guard for a known
old-miniprot long-sequence defect. It does not claim to have reproduced the
reporter's private 20-Gb data or host state.

## Evidence boundary

The issue log establishes the following facts:

- The target contained 20,029,007,188 bases in 657 sequences.
- LiftOn was invoked with 40 threads and copy search enabled.
- minimap2 and miniprot each returned Python subprocess code `-11`, which is
  POSIX signal 11 (`SIGSEGV`).
- miniprot's last reported milestone was `collected syncmers`.
- Experimental GTF conversion was followed by 101 input-validation errors.

Those observations do not establish why either native process received
`SIGSEGV`. In particular, `-11` is not proof of an OOM kill. Confirming the
original cause still requires the input manifest and native-tool versions,
plus scheduler, cgroup, kernel, and process-memory evidence from the
reporter's run.

The 101 annotation findings are a separate correctness risk: the converted
GFF3 and `stats/gff3_input_validation.txt` must be reviewed before accepting a
successful rerun. The log provides no evidence that those findings caused the
two independent native `SIGSEGV` terminations.

The exact v1.0.11 release remains tag `v1.0.11` at commit `c623f0b`. Its sealed
scientific evaluation is not reinterpreted or amended by this release. The
post-tag trans-spliced-copy repair at `c8daa44` is included in v1.0.12 and is
documented as post-v1.0.11 work, never as behavior measured in the exact-tag
campaign.

## Controlled mechanism reproduction

Core dumps were disabled for all controls. On the repository's 50,818,468-bp
chromosome-22 fixture:

| Native command | Unconstrained result | Constrained result |
|---|---:|---:|
| miniprot 0.13-r248 | exit 0; 0.598 GB reported peak RSS | `prlimit --as=943718400`: `-11`, after `collected syncmers` |
| minimap2 2.28-r1209 | exit 0; 0.390 GB reported peak RSS | `prlimit --as=419430400`: `-11` |

Lower address-space ceilings caused allocator-explicit aborts in some
miniprot controls, showing that resource pressure need not have one invariant
failure signature. These experiments reproduce a plausible mechanism and the
issue-specific miniprot milestone/signature combination. They do not recover
the reporter's causal history.

## Public 22-Gb scale surrogate

Because the issue inputs are private, scale behavior was measured on the NCBI
[Ptaeda2.0 loblolly-pine assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000404065.3/)
(`GCA_000404065.3`):

- 22,103,635,615 bases in 1,760,464 sequences; maximum sequence 1,058,009 bp;
  N50 107,038 bp; 1,571,989,997 `N` bases.
- Compressed NCBI-file MD5: `0034bdef1a866ef20b3d981bcebee6c0`.
- Decompressed FASTA SHA-256:
  `f81866cd46cb2d0b9517f22c18c7342233778916ee964c8d7873c85eacdb9f0f`.
- DNA query (`test/test.fa`) SHA-256:
  `319b94d37578083c9a4d6d555e838431ae4e9c8758ceba4200a8bdff573e3512`.
- Protein query (`test/test_prot.fa`) SHA-256:
  `c156cc04b01ec942483a2de01ace2368d4ff268a9edc4abf756905fd98ec9a98`.

The host exposed 1,081,490,325,504 bytes of RAM and no swap. A watchdog
sampled process-group RSS every second, capped the group at 700 GiB, reserved
at least 200 GiB for the host, imposed a 24-hour timeout, and disabled core
dumps. All native commands used 40 threads and LiftOn's corresponding native
options.

Native-tool provenance was frozen as follows:

- installed miniprot 0.13-r248 binary SHA-256
  `f7da26eb8843fc33fb2f326b0ce9259d251247dc971100792ba3e3db28cbd38d`;
- installed minimap2 2.28-r1209 binary SHA-256
  `56e910107685c2448c5b673e094098bfdf70472bfa46f0b38c617edb8438e0dd`;
- source-built miniprot 0.18 commit
  `671db243f964a68bd724af11cd9964d840f29c43`, binary SHA-256
  `1331266b654957ea2efb5509c2007f582444040fd71eafd3211deb26a7be1bc2`;
- source-built minimap2 2.31 commit
  `3c28777e7e2dcc90f825de1b9f17a89cca7d4452`, binary SHA-256
  `7e6e748bfba4cb72b73d168fa4e3ba717d0561a84e9ef15e5cc3be0322542ecb`.

| Program | Version | Exit | Wall seconds | Peak group RSS |
|---|---|---:|---:|---:|
| miniprot | 0.18-r281 | 0 | 830.862 | 119,382,298,624 B (111.18 GiB) |
| miniprot | 0.13-r248 | 0 | 549.154 | 119,520,882,688 B (111.31 GiB) |
| minimap2 | 2.31-r1302 | 0 | 324.554 | 59,006,746,624 B (54.95 GiB) |
| minimap2 | 2.28-r1209 | 0 | 345.890 | 56,678,670,336 B (52.79 GiB) |

The installed miniprot 0.13/minimap2 2.28 pair also completed when started
concurrently. It took 629.022 seconds and peaked at 161,202,130,944 bytes
(150.13 GiB). Sequential execution bounds the pair's combined peak by the
larger solo maximum, 111.31 GiB: 38.82 GiB (25.9%) below the observed
concurrent peak. The one-second trace SHA-256 is
`3600a80534a7793ddcf8963cf336c7564e7405475c64db4d6117a423b0fd8e6c`.
At the peak sample (251.24 seconds), the trace assigned 102,436,405,248 bytes
to the miniprot process group and 58,763,628,544 bytes to minimap2; the wrapper
shell accounted for the small remainder.

This surrogate is much more fragmented than the reported target and uses one
small query to isolate target indexing. It does not exercise the reporter's
74,405 proteins, subsequent annotation processing, or unknown host limits.
Successful completion on a 1-TiB host therefore does not explain the private
failure. It does demonstrate that both native indexes are individually large
at this genome scale and that overlap can materially increase peak memory.

## v1.0.12 contracts

- The default remains concurrent at or below 4,000,000,000 target bases and
  becomes sequential above that boundary. This reuses LiftOn's existing
  minimap2 split-index boundary; it is not a per-host RAM guarantee.
- `--serial-aligners` and `--parallel-aligners` are explicit, mutually
  exclusive overrides. A forced large-target parallel run warns about peak
  overlap.
- Liftoff creates at most one Python worker per real alignment task and divides
  the requested thread budget among workers. A one-task, 40-thread mapping no
  longer creates 40 mostly idle workers.
- The run manifest records target size, longest sequence, resolved schedule,
  exact native commands, return/signal facts, last detected stage, and no more
  than 64 KiB of stderr per execution.
- Diagnostics report facts and never infer OOM from a signal.
- A target sequence at least `2^31` bases long is rejected with a parseable
  miniprot version older than 0.14. An unparseable version warns. The threshold
  follows miniprot's [upstream 32-to-64-bit sequence-length
  correction](https://github.com/lh3/miniprot/commit/ec4fdba0f01485e3fce9874a9ac1c3e40d409543).
- Scheduling changes do not alter thresholds, candidate scoring, rescue, or
  best-outcome selection.

## Local qualification gates

| Gate | Result |
|---|---|
| Focused issue/regression matrix | **112 passed** |
| Fatal flake8 rules over repository | **passed** |
| Sphinx HTML build | **passed**; 78 pre-existing warnings, no new warning class |
| `make test-fast LIFTON_PY=...` | **30 passed** |
| `make benchmark-gate LIFTON_PY=...` | **PASS**; exact quality; protein identity 0.99425; coding and total completeness 0.99764; wall time 44% below baseline |
| Complete `pytest tests/ -q` (isolated worktree, `PYTHONHASHSEED=0`) | **1,856 passed, 2 skipped** (2026-09-11) |
| Real chr22 default/serial byte identity (minimap2 2.28, miniprot 0.13, `-copies -sc 0.95`) | **byte-identical**, md5 `160e05449efdb8bccef06fa12bb9b8b4`; 2 min 36 s each; `gff3-validate` VALID, 0 errors, 59 warnings |
| Wheel + sdist build and clean-wheel smoke | **passed**: `python -m build`; `pip install --no-cache-dir` of the wheel into a fresh venv; `lifton --version`, `lifton -h`, `gff3-validate -h`, and `import lifton.tool_execution` all succeed |

## Scope added on 2026-09-11

After the #71 work was gated and committed (`307abc6`), v1.0.12 took on the
improvement program in `notes/lifton_v1.0.12_improvement_analysis.md`. The #71
contract above, that scheduling changes alter no thresholds or scoring, still
describes the #71 changes. It no longer describes the release: v1.0.12 now also
changes rescue output by default.

| Commit | Change | Output | Evidence |
|---|---|---|---|
| `e282103` | gene-level, primary-assembly, GeneID recall; recall-gap diagnostic | benchmark tooling only | unit tests; re-derived the analysis numbers |
| `16f67d2` | forked workers can no longer abort the parent's output transaction | none (latent fault) | fork regression test, fails on old code |
| `79a212a` | protein-coverage rescue sub-pass + isoform-aware rescue, default on | **adds genes and transcripts** at distance | strict A/B 13/13 + 13/13: 0 lost, 0 duplicates, 0 regressions, validity unchanged |
| `d4c1f02` | `--locus-pipeline` by default when `-t > 1` | byte-identical | drosophila and dog → cat whole genomes byte-identical; wall −33 % |
| `764e1eb` | intermediate files back under `lifton_output/` | paths only | regression tests, fail on old code |
| `8816ef3` | faster Liftoff SAM parsing and GFF3 writing | byte-identical | 120,000-case equivalence; identical Liftoff bodies on two whole genomes |
| `804f01e` | forked parallel Liftoff lift loop, default when `-t > 1` | byte-identical | identical on drosophila and dog → cat (4 arms each); aligner phase −32 to −33 % |

Remaining release gates for the combined tree:

| Gate | Result |
|---|---|
| Complete `pytest tests/` on the final tree | **1,924 passed, 2 skipped** at `804f01e` (isolated worktree) |
| Real chr22 example with the new defaults (`-copies -sc 0.95`) | exit 0; `-t 1` (2 min 32 s) and `-t 8` (1 min 23 s) byte-identical, and identical to the `307abc6` output (same-species: no rescue candidates); `gff3-validate` VALID, 0 errors, 59 warnings; intermediates in `lifton_output/liftoff/` and `lifton_output/miniprot/` |
| Wheel + sdist build and clean-wheel smoke | **passed** at `804f01e`: `--no-cache-dir` install into a fresh venv; `lifton --version` v1.0.12; new opt-outs in `-h`; promoted defaults active; `gff3-validate` present |
| `devel` CI on Python 3.10/3.11/3.12 | pending push |
| Tag / GitHub Release / PyPI | pending sign-off |

## Scope added on 2026-09-12 (second round)

The first round's results made a second measurement possible: replaying the
gates LiftOn now *ships* against its own whole-genome output. It said the
distant-recall well is nearly dry — 79–91 % of what is still missed is a locus
another gene already holds — and that the weak point had moved to model
quality. The round that followed is recorded in §5.5–5.7 of
`notes/lifton_v1.0.12_improvement_analysis.md`.

| Commit | Change | Output | Evidence |
|---|---|---|---|
| `adfc9a0` | track the Figure-4 diagnostic `cross_locus_rescue.py` cites | benchmark tooling only | inventory |
| `f43e044` | recall-gap diagnostic replays the shipped gates; rescued-model ORF validity; co-ortholog count | benchmark tooling only | re-derived by hand first, then by the tool |
| `d3fa4da` `cff6ec0` `8b91f3a` | terminal-stop completion for miniprot-derived models, default on | **changes miniprot-derived CDS** by exactly one stop codon | ladder A/B 8/8: 0 changed other than a three-base terminal extension, 0 regressions |
| `ec71375` | flat annotations lift (#37); `-dir/--intermediate-dir` (#14) | corrective / paths only | regression tests that fail on old code |
| `e98ef13` | cross-locus replacement keeps its isoforms, gates on coverage, keeps reference transcript ids | opt-in pass only | unit + end-to-end tests; A/B pending |
| `2a8fee0` | isoform prefetch stops re-fetching a row it already held | byte-identical | two subsets, `-t 1` and `-t 8`, identical md5 |
| `843a73e` `2a91017` | analysis note, changelogs, CLAUDE.md | docs | — |

### Dependabot: three open pyo3 advisories, assessed and deferred

`lifton/gffbase/_rust/Cargo.lock` pins pyo3 0.22.6, and three advisories are
open against it (one high). They are **not release blockers**, for three
reasons, and the fix is not a release-time change:

- the crate uses none of the three affected APIs — no `PyString::from_object`,
  no `nth`/`nth_back` on `PyList`/`PyTuple`, no `PyCFunction::new_closure`;
- the Rust extension is not built here and is not on the default path — a
  missing `_native` falls back to the pure-Python parser, which is what the
  suite exercises;
- two of the three are first patched in pyo3 **0.29.0**, seven majors ahead of
  the pin, and the toolchain available for validation is Rust 1.69.

Recorded as follow-up work with its own migration and a rebuilt extension, not
as something to attempt while cutting a release.

## Release and deployment surfaces

| Surface | State |
|---|---|
| Git commit | pending |
| `devel` CI on Python 3.10/3.11/3.12 | pending |
| `main` CI and Sphinx Pages | pending |
| Git tag / GitHub Release | pending |
| TestPyPI / PyPI | pending |
| Issue #71 follow-up | pending reporter confirmation after release |
