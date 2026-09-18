# LiftOn v1.0.13 packaging release — issue #78

Approved 2026-09-16. Released baseline: v1.0.12, 6c86d1b.
Upstream synchronized to 8ebf43b7d7170dc32f0f66d0ad698c4b43e32533.
Implementation: `/tmp/lifton-issue78-release`, `release-1.0.13-packaging`.
Earlier correctness work remains at `/tmp/lifton-v1.0.13-build`, 9bdf8a1;
its branch name is historical and does not allocate the packaging release.
Resume that work as the subsequent correctness milestone after this hotfix.

## Required behavior

- Standard pip install excludes mappy. `lifton[native]` adds the optional binding;
  activation still requires `--native` and `LIFTON_NATIVE_LIFTOFF_ALIGN=1`.
- Missing binding preserves the existing subprocess fallback and clear warning.
- Pip does not supply minimap2/miniprot executables. Fresh standard lifts require
  both on PATH; evaluation and valid precomputed alignments retain their bypasses.
- Bioconda's noarch Python package supplies prebuilt mappy and both aligners.
- Ordinary annotation behavior remains v1.0.12. No algorithm changes, broad
  vendored refactors or speculative dependency upgrades belong to this release.
- No publishing, pushing, deployed documentation or external issue/PR messages.

## Execution ledger

| Step | Status | Gate |
|---|---|---|
| Synchronize and isolate | Complete | Upstream/release/issue/PR checked; clean separate worktree. |
| Optional mappy and fallback | Complete | 89 focused tests pass; fallback verified with mappy absent; warnings actionable. |
| External tools and documentation | Complete | Actionable preflight and corrected pip/source/macOS/Seqera commands. |
| Seqera/compiler-free installs | Complete | Python 3.10–3.12 and 3.14 compiler-free Singularity containers pass all 8 checks. |
| Bioconda recipe | Complete | Update PR #66594 recipe locally; cigar/cap removed, dependencies aligned, bioconda-utils lint passes. |
| Frozen qualification | Complete | Full suite passes (2,043 passed × 3 interpreters, 0 failed), 24 cells byte-identical. |
| Release packet | Complete | Exact artifacts/evidence, recipe overlay, draft responses and explicit public gates. |

Intermediate validation: regression tests first reproduced three intended failures.
Focused checks: 89 passed in 1.96s. Fast/native matrix: 30 passed (all 24 cells).
Fatal flake8 and diffcheck pass. Draft wheel fresh no-cache install succeeds with
mappy absent. Constructed two-strand native lifts recover exact independent proteins
and all gene/mRNA/exon/CDS rows, ordinary and stream/inmemory byte-identical.
Draft evidence: `/tmp/lifton-issue78-evidence/pilot-native/evidence.json`.
These checks precede the frozen artifact/full-suite qualification.

CI changes: reusable packaging workflow builds the exact wheel/sdist, tests each
in compiler-free Python 3.10–3.12/3.14 containers and runs fresh native lifts.
Publishing depends on completion of that workflow. Core tests retain the complete
suite and all 24 configurations; optional binding tests have a separate matrix.
No remote workflow has been triggered and no publication is authorized.

## Evidence policy and release ordering

Fresh environments outside source, no cache for package installs, exact imports,
input/source/tool/dependency/artifact hashes and commands. Full suite once per
frozen candidate on Python 3.10–3.12; binding tests separately with native extra.
Python 3.14 receives installation/smoke checks, not a new scientific-support claim.
Native comparisons hold aligner versions and options fixed; successful installation
alone does not establish a usable annotation container. Small constructed inputs
have independent sequence/coordinate truth; real GFF3 control checks parity.

Local Docker/Podman daemons are unavailable. Use Singularity/isolated Conda locally
and compiler-free container CI. Local recipe uses a verified candidate-sdist source
overlay; final public PyPI URL/checksum is verified after authorized publication.
Public Seqera request and refreshed Bioconda CI remain separate gates, never inferred
from older builds. At most 32 aggregate native threads and 256 GiB host reserve.

## Primary sources

- [Issue #78](https://github.com/Kuanhao-Chao/LiftOn/issues/78), including miniprot comment.
- [Reported Wave build](https://wave.seqera.io/view/builds/bd-6ea91f07506bfbd7_2):
  Python 3.14.7, mappy 2.31 sdist, missing gcc; no need to compile binding for default path.
- [Bioconda PR #66594](https://github.com/bioconda/bioconda-recipes/pull/66594):
  still v1.0.9; review label already present; historical CI passed.
- [Bioconda guidelines](https://bioconda.github.io/contributor/guidelines.html)
  and [local testing](https://bioconda.github.io/contributor/building-locally.html).
