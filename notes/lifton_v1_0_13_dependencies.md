# LiftOn v1.0.13 dependency audit

Checked 2026-09-14. This is an audit, not an upgrade proposal. Scientific
comparisons retain the tested environment and record content hashes. No runtime
requirements have been changed by this audit.

| Package | Installed | Current upstream | Upstream Python requirement |
|---|---|---|---|
| lifton | 1.0.10 | [1.0.12](https://pypi.org/project/lifton/1.0.12/) | >=3.10 |
| numpy | 2.4.4 | [2.5.3](https://pypi.org/project/numpy/2.5.3/) | >=3.12 |
| biopython | 1.87 | [1.88](https://pypi.org/project/biopython/1.88/) | >=3.10 |
| parasail | 1.3.4 | [1.3.4](https://pypi.org/project/parasail/1.3.4/) | unspecified |
| intervaltree | 3.2.1 | [3.2.1](https://pypi.org/project/intervaltree/3.2.1/) | !=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*,!=3.4.*,>=2.7.18 |
| interlap | 0.2.7 | [0.2.7](https://pypi.org/project/interlap/0.2.7/) | unspecified |
| networkx | 3.6.1 | [3.6.1](https://pypi.org/project/networkx/3.6.1/) | !=3.14.1,>=3.11 |
| pyfaidx | 0.9.0.4 | [0.9.0.4](https://pypi.org/project/pyfaidx/0.9.0.4/) | >=3.7 |
| pysam | 0.24.0 | [0.24.1](https://pypi.org/project/pysam/0.24.1/) | >=3.9 |
| gffutils | 0.14 | [0.14](https://pypi.org/project/gffutils/0.14/) | >=3.8 |
| ujson | 5.12.0 | [6.0.0](https://pypi.org/project/ujson/6.0.0/) | >=3.10 |
| duckdb | 1.5.2 | [1.5.5](https://pypi.org/project/duckdb/1.5.5/) | >=3.10.0 |
| pyarrow | 24.0.0 | [25.0.1](https://pypi.org/project/pyarrow/25.0.1/) | >=3.10 |
| mappy | 2.30 | [2.31](https://pypi.org/project/mappy/2.31/) | unspecified |
| pytest | 9.0.3 | [9.1.1](https://pypi.org/project/pytest/9.1.1/) | >=3.10 |
| hypothesis | 6.155.2 | [6.168.0](https://pypi.org/project/hypothesis/6.168.0/) | >=3.10 |
| coverage | absent | [7.16.1](https://pypi.org/project/coverage/7.16.1/) | >=3.10 |
| flake8 | 7.3.0 | [7.3.0](https://pypi.org/project/flake8/7.3.0/) | >=3.9 |
| build | 1.5.0 | [1.6.1](https://pypi.org/project/build/1.6.1/) | >=3.10 |
| setuptools | 82.0.1 | [84.0.0](https://pypi.org/project/setuptools/84.0.0/) | >=3.10 |
| pip | 26.0.1 | [26.2.1](https://pypi.org/project/pip/26.2.1/) | >=3.10 |

The shared environment advertises installed LiftOn 1.0.10, whereas importing
from the release checkout yields v1.0.12. Qualification must record actual
module path/version and source commit; installed distribution metadata alone
does not identify executed source. Coverage is absent from the current shared
environment and will be installed in the isolated qualification environment.

Native tools in the current scientific environment:

| Tool | Installed | Current upstream tag | SHA-256 of installed executable |
|---|---|---|---|
| minimap2 | 2.28-r1209 | [v2.31](https://github.com/lh3/minimap2/releases/tag/v2.31) | `56e910107685c2448c5b673e094098bfdf70472bfa46f0b38c617edb8438e0dd` |
| miniprot | 0.13-r248 | [v0.18](https://github.com/lh3/miniprot/releases/tag/v0.18) | `f7da26eb8843fc33fb2f326b0ce9259d251247dc971100792ba3e3db28cbd38d` |
| gffread | 0.12.7 | [not checked](https://github.com/gpertea/gffread) | `824cf258d9e856261253aa149c1617517e64d7ab92d112fe275853783a866a9d` |
| gt | /ccb/sw/bin/gt (GenomeTools) 1.6.1 | [not checked](https://genometools.org/) | `e6f553e445a1cb5fd7078139c39fe0738a5cfd7058a148edf41a65e784078fbd` |

## Decisions and qualification implications

- Keep native minimap2 2.28-r1209 and miniprot 0.13-r248 fixed during the
  correctness comparison. Newer native releases can change placements and need
  a separate paired experiment before adoption.
- Keep DuckDB 1.5.3/1.5.4 exclusions. The repository records a large streamed
  GEOMETRY append failure (GH #56); a newer upstream release does not by itself
  prove that LiftOn-specific regression fixed. Test 1.5.5 separately before any
  exclusion/fallback changes.
- Latest NumPy/PyArrow or other packages can raise the Python floor. Use
  interpreter-compatible pinned qualification environments; do not force a
  Python 3.11 lock onto Python 3.10. Record each resolved environment separately.
- Build-system requirements in pyproject.toml use version ranges, despite a
  comment calling the backend pinned. Record actual build versions and artifact
  hashes; do not describe those ranges as a reproducible lock.
- Preserve unchanged requirements unless a demonstrated correctness or install
  failure requires an isolated dependency fix. Plain package installation and
  a scientific benchmark environment are separate qualification checks.

Machine-readable local and upstream metadata are retained in the execution
workspace; final campaign receipts will seal dependency file hashes, native
binaries, interpreter, sources and inputs. No claim of vulnerability scanning
is made by this version audit.
