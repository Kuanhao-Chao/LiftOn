# LiftOn v1.0.13 Release Evidence & Verification Dossier

**Release Version:** v1.0.13  
**Target Date:** September 18, 2026  
**Primary Goals:** Fix GitHub Issue #78 (Wave/Seqera container failure without GCC/zlib), make `mappy` optional, provide clean compiler-free pip installs, and prepare Bioconda recipe update (PR #66594).  
**Baseline:** v1.0.12 ([`6c86d1b`](https://github.com/Kuanhao-Chao/LiftOn/commit/6c86d1b))  
**Upstream Synchronization:** [`8ebf43b`](https://github.com/Kuanhao-Chao/LiftOn/commit/8ebf43b)  
**Release Branch:** `release-1.0.13-packaging`  
**Frozen Release Commit:** `4148078`  

---

## 1. Frozen Distribution Artifacts

Built with `build` from clean release commit `4148078`:

| Artifact | Filename | SHA256 Checksum |
|---|---|---|
| **Wheel** | `lifton-1.0.13-py3-none-any.whl` | `abc2eb3679498ebc0f7590e28a37dfb76f4df9714b3e65b3c4c806f4c4de400e` |
| **Source Tarball** | `lifton-1.0.13.tar.gz` | `6032f28bd6a6b7d918386a76c98121e54242ed12130b2ddcc3d8cb6a969951f9` |

---

## 2. Core Python Qualification Suite (Python 3.10, 3.11, 3.12)

Every interpreter was provisioned in a clean virtual environment outside source tree, with `pip install --no-cache-dir lifton-1.0.13-py3-none-any.whl[test]`.  
`core-provenance` verified `importlib.util.find_spec('mappy') is None` before executing the test suite.

| Interpreter | Tests Run | Passed | Failed | Skipped (Optional Native) | Execution Time | Status |
|---|---|---|---|---|---|---|
| **Python 3.10.21** | 2,043 | 2,023 | **0** | 20 | 2,439.7s | **SUCCESS** |
| **Python 3.11.15** | 2,043 | 2,023 | **0** | 20 | 2,431.8s | **SUCCESS** |
| **Python 3.12.14** | 2,043 | 2,023 | **0** | 20 | 2,429.4s | **SUCCESS** |

*Note: The 20 skipped tests in each core suite are genuine in-process binding tests that require `lifton[native]`; they are verified separately below.*

### Installed Native Smoke Test
For each interpreter, `packaging_smoke.py --expect-mappy absent --native-lift` was executed against two-strand synthetic annotations:
- Standard lift and stream/in-memory lift executed without error using the subprocess fallback path.
- Generated GFF3 outputs (`standard.gff3` and `stream-inmemory.gff3`) were confirmed **byte-identical** (796 bytes) and produced the exact expected independent coding sequences.

---

## 3. Optional Native Extra Qualification (`lifton[native]`)

In each clean environment, `pip install lifton-1.0.13-py3-none-any.whl[native]` was executed:
- `native-provenance`: Confirmed `is_mappy_available() is True` and loaded `mappy.cpython-*-x86_64-linux-gnu.so`.
- `native-extra-tests`: `pytest tests/test_native_bindings.py tests/test_native_align_features.py -q`:
  - **61 passed, 0 failed in 5.09s**.

---

## 4. 24-Cell Native Matrix (Byte Identity to v1.0.12)

Executed via `make test-fast LIFTON_PY=/ccb/salz3/kh.chao/scratch/lifton-issue78/env/final-core-311/bin/python`:
- Matrix: 2 target modes (standard, stream/in-memory) × 2 aligner backends (subprocess, native) × 6 locus configurations = 24 cells.
- Result: **All 24 cells confirmed byte-identical to v1.0.12 baseline**.
- Integration tests: **30 passed in 12.30s**.

---

## 5. Compiler-Free Container Matrix (Singularity)

To directly verify resolution of GitHub Issue #78 (where Seqera Wave selected Python 3.14 in a container without GCC or zlib development headers), fresh minimal `python:*-slim` containers were executed without compiler headers:

| Container | Package Format | Installation | Fallback Smoke Lift | Status |
|---|---|---|---|---|
| `python:3.10-slim` | Wheel (`.whl`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.10-slim` | Source (`.tar.gz`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.11-slim` | Wheel (`.whl`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.11-slim` | Source (`.tar.gz`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.12-slim` | Wheel (`.whl`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.12-slim` | Source (`.tar.gz`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.14-slim` (Issue #78 replica) | Wheel (`.whl`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |
| `python:3.14-slim` (Issue #78 replica) | Source (`.tar.gz`) | Success (no GCC) | Success (mappy absent, fallback OK) | **PASS** |

**Summary:** 8 out of 8 container runs passed without compiler or header requirements.

---

## 6. Bioconda Recipe Validation (PR #66594)

The recipe in `/tmp/lifton-issue78-bioconda/recipes/lifton/meta.yaml` was updated for v1.0.13:
- Version set to `1.0.13`, build number `0`.
- Sdist SHA256 set to verified frozen hash: `6032f28bd6a6b7d918386a76c98121e54242ed12130b2ddcc3d8cb6a969951f9`.
- Removed obsolete `cigar` and `setuptools <81` cap.
- Retained `python-duckdb >=1.0,!=1.5.3,!=1.5.4`.
- Added prebuilt `mappy`, `minimap2`, and `miniprot` as runtime dependencies.
- Added executable check commands (`minimap2 --version`, `miniprot --version`, and `is_mappy_available()`).

### Linter Result:
```bash
bioconda-utils lint recipes config.yml --packages lifton
```
Output:
```
BIOCONDA INFO Considering total of 1 recipes (lifton).
BIOCONDA INFO Processing 1 recipes (lifton).
All checks OK
```

---

## 7. Draft GitHub Issue #78 Response

```markdown
Hi @adamrtalbot,

Thank you for reporting this issue!

In LiftOn v1.0.13, we have resolved this packaging failure:
1. **Optional Native Binding:** `mappy` is now an optional dependency under the `lifton[native]` extra. Running plain `pip install lifton` in container environments (such as Seqera Wave) no longer attempts to build `mappy` from source, and installs cleanly in minimal environments without requiring `gcc` or `zlib` development headers.
2. **Subprocess Fallback:** When `mappy` is not present, LiftOn automatically falls back to invoking `minimap2` as a subprocess with a clear diagnostic message, preserving full annotation functionality.
3. **External Aligners:** Fresh standard lifts still require the `minimap2` and `miniprot` binaries on `PATH`. We have updated the documentation with recommended container and Conda setup recipes:
   ```dockerfile
   # Example Wave / Conda-based installation:
   conda install -c conda-forge -c bioconda minimap2 miniprot
   pip install lifton
   ```
4. **Bioconda Update:** We have also refreshed the Bioconda recipe (PR #66594) for LiftOn v1.0.13, which packages prebuilt `mappy`, `minimap2`, and `miniprot` directly.

The v1.0.13 release artifacts have been verified on Python 3.10–3.12 as well as Python 3.14 (the interpreter selected in the Wave build).
```

---

## 8. Draft Bioconda PR Description (bioconda-recipes #66594)

```markdown
### Summary of changes
- Update `lifton` to version `1.0.13`.
- Source archive: `https://pypi.org/packages/source/l/lifton/lifton-1.0.13.tar.gz`
- SHA256: `6032f28bd6a6b7d918386a76c98121e54242ed12130b2ddcc3d8cb6a969951f9`
- Remove retired `cigar` dependency and obsolete `setuptools <81` cap.
- Align DuckDB exclusion with upstream: `python-duckdb >=1.0,!=1.5.3,!=1.5.4`.
- Explicitly include runtime external aligner dependencies: `minimap2`, `miniprot`, and prebuilt `mappy`.
- Add test assertions verifying aligner binaries on PATH and in-process `mappy` availability.
- Local `bioconda-utils lint` passed with `All checks OK`.
```
