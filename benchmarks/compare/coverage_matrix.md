# LiftOn website-example coverage matrix

Audit of the 8 examples on https://khchao.com/LiftOn/ against the subset registry (`benchmarks.json`), the full-genome FTP registry (`datasets.json`), and on-disk inputs.

| # | website example | class | subset id | in benchmarks.json | full id | in datasets.json | ref inputs on disk | subset built (S/L/M) |
|---|---|---|---|---|---|---|---|---|
| 1 | Human GRCh38 → T2T-CHM13 | same-species | `human_mane` | ✅ | `human` | ✅ | ✅ | ✅✅✅ |
| 2 | Mouse GRCm39 → NOD/SCID | same-species | `mouse` | ✅ | `mouse` | ✅ | ✅ | ✅✅✅ |
| 3 | Honey bee HAv3.1 → ASM1932182v1 | same-species | `bee` | ✅ | `bee` | ✅ | ✅ | ✅✅✅ |
| 4 | Arabidopsis TAIR10 → ASM2311539v1 | same-species | `arabidopsis` | ✅ | `arabidopsis` | ✅ | ✅ | ✅✅✅ |
| 5 | Rice IRGSP → ASM3414082v1 | same-species | `rice` | ✅ | `rice` | ✅ | ✅ | ✅✅✅ |
| 6 | Human GRCh38 → Chimpanzee | close cross-species | `human_to_chimp` | ✅ | `human_to_chimp` | ✅ | ✅ | ✅✅✅ |
| 7 | D. melanogaster → D. erecta | distant cross-species | `drosophila` | ✅ | `drosophila` | ✅ | ✅ | ✅✅✅ |
| 8 | Mouse GRCm39 → Rat mRatBN7.2 | distant cross-species | `mouse_to_rat` | ✅ | `mouse_to_rat` | ✅ | ✅ | ✅✅✅ |

**All 8 website examples covered by a subset benchmark with on-disk inputs:** YES.

## Full-genome FTP registry (`datasets.json`) membership

| website full id | present in datasets.json |
|---|---|
| `human` | ✅ |
| `mouse` | ✅ |
| `bee` | ✅ |
| `arabidopsis` | ✅ |
| `rice` | ✅ |
| `human_to_chimp` | ✅ |
| `drosophila` | ✅ |
| `mouse_to_rat` | ✅ |

Missing from `datasets.json`: none — registry complete.

