# Release readiness — the 27 unpushed commits on `devel`

Written 2026-07-26. **Nothing has been pushed or tagged, and `__version__` is
untouched at `v1.0.10.1`.** This is an audit and a recommendation, not a release.

## What is on the branch

27 commits ahead of `origin/devel`, in four batches:

| batch | commits | character |
|---|---|---|
| v1.0.10.1 issue fixes | `23b4dc1`…`0a0d547` | 8 GitHub issues |
| full-repository audit | `c8ef69b`…`c1e44e9` | 8 correctness fixes, 3 perf wins, 1 self-caught NO-GO |
| CDS attribute parity | `2085061`…`c54c396` | 1 output change, 4 robustness fixes, 2 perf items |
| `start_lost` diagnostic | `7ce9d45` | NO-GO, documented |

Suite **1668 passed / 2 skipped**. The 24-cell matrix is green with no golden
edit throughout, and the drosophila anchor reproduces both of its md5s
(`3f908d02…` under `LIFTON_NO_CDS_ATTR_CARRY=1`, `396fd60b…` by default).

## Do the committed benchmark figures still hold?

**Accuracy and completeness: YES, provably — verified, not assumed.**

Two independent lines of evidence:

1. On the full dog→cat genome, **columns 1-8 are byte-identical across all
   1,815,199 rows**. For a CDS those columns *are* its coordinates, strand and
   phase, so no lifted model moved.
2. `benchmarks/compare/evaluator.py` scores by re-aligning sequence fetched from
   those coordinates and computes `protein_identity` itself. It never reads
   `Dbxref`, `product`, `protein_id`, `gbkey` or `locus_tag` — the attributes
   this batch added.

Run directly to confirm rather than infer: the neutral evaluator over the
drosophila off/on arms produced an **identical `transcripts.tsv` and
`feature_types.tsv`**, with `summary.json` differing only in the arm's own name
and output path. Every metric matched.

**Validity: unchanged.** `gff3-validate` reports 0 errors / 176 warnings on the
full dog→cat output in both arms, and 0 / 66 on drosophila in both arms.

**Performance: this is the one panel a release should refresh.** Output grows
12–43 % (drosophila 29.4 → 38.3 MB; dog→cat 398 → 524 MB), which moves write
time and downstream file sizes, and Step-7 wall improved measurably on mammalian
data (−5.4 %, faster in 6 of 6 alternated repeats). Peak RSS moved +2.5 %.

`benchmarks/compare/fourway_results.json` is frozen behind a SHA guard in
`benchmarks/inventory_rules.json`, so any refresh is a deliberate, test-enforced
act — as it should be.

## What a release would need

1. **A version decision.** The batch contains output-affecting changes (CDS
   attribute parity; copy-gene CDS identifiers), so it is not a patch release
   over `v1.0.10.1`. **v1.0.11** fits: outputs differ from v1.0.10.1 in column 9
   and in copy-gene identifiers, while coordinates and scores do not.
2. **A perf-panel refresh only.** Accuracy, completeness and validity are
   unchanged and need no re-baseline; wall/memory/file-size figures do.
3. **A note for downstream consumers.** Output files are 12–43 % larger. Anyone
   parsing CDS rows now sees the reference's descriptive attributes where before
   they saw only `ID`/`Parent`; `LIFTON_NO_CDS_ATTR_CARRY=1` restores the old
   shape for reproduction of earlier results.

## Deliberately NOT done here

No `__version__` bump, no tag, no push, no `fourway_results.json` edit. All of
those need sign-off.
