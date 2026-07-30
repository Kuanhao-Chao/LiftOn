# Release readiness — v1.0.10

Written 2026-07-30. **Nothing has been pushed and no tag has been moved.** This
records what the release contains, what was verified, and the exact sequence that
publishes it.

## Why v1.0.10 and not v1.0.11

`v1.0.10` was tagged at `3739dfc` and never released. There is no GitHub Release
for it (the latest is v1.0.9), nothing on PyPI (the latest is `1.0.9`), and the
tag never reached `main` — v1.0.8 and v1.0.9 are both reachable from `main` and
`devel`. Meanwhile `__version__` drifted to `v1.0.10.1`, a version prepared for a
GitHub-issue batch and likewise never published.

So there were two version numbers users could read about and neither could be
installed. Rather than leave both stranded and number the next release v1.0.11,
`v1.0.10` now means the whole body of work: the originally staged v1.0.10, the
issue batch, the full-repository correctness audit, and the CDS-attribute parity
fix. **v1.0.9 is therefore the baseline for every claim in `CHANGELOG.md`.**

Moving the tag is safe precisely because nothing depends on where it points: no
Release, no PyPI artifact, no `main` ancestry.

## What is on the branch

67 commits past the `v1.0.10` tag, in five batches:

| batch | commits | character |
|---|---|---|
| v1.0.10.1 issue fixes | `23b4dc1`…`0a0d547` | 8 GitHub issues |
| release qualification / robustness | `d503069`…`e7a2aeb` | manifests, transactional output, caches, bounded scheduling |
| full-repository audit | `c8ef69b`…`c1e44e9` | 8 correctness fixes, 3 perf wins, 1 self-caught NO-GO |
| CDS attribute parity | `2085061`…`c54c396` | 1 output change, 4 robustness fixes, 2 perf items |
| streaming validator / SQL / release prep | `7ce9d45`…HEAD | 21× validator memory, −9.6% SQL, this release |

## Do the committed benchmark figures still hold?

**Accuracy and completeness: yes, provably — verified, not assumed.**

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

**Performance: this is the one panel a benchmark refresh should update.** Output
grows 12–43 % (drosophila 29.4 → 38.3 MB; dog→cat 398 → 524 MB), which moves
write time and downstream file sizes, and Step-7 wall improved measurably on
mammalian data (−5.4 %, faster in 6 of 6 alternated repeats). Peak RSS moved
+2.5 %.

`benchmarks/compare/fourway_results.json` is frozen behind a SHA guard in
`benchmarks/inventory_rules.json`, so any refresh is a deliberate, test-enforced
act — as it should be.

## Verification performed for this release

See the release-preparation session for the full log. Summary of the gates:

- `flake8 . --select=E9,F63,F7,F82` — clean
- Full `pytest tests/` suite
- The 24-cell byte-identity matrix, green with no golden edit
- `make test-fault` — injected write, validation, cache and subprocess failures
- `python -m build`, `twine check`, and the wheel installed into a clean venv
  (`lifton -V` prints `v1.0.10`)
- The chr22 shell example against real minimap2/miniprot, then `gff3-validate`
  on its output
- The drosophila byte anchors, with cached `-L`/`-M` at `-t 1`

## Publishing (not done — requires sign-off)

```bash
# 1. Push both branches. main is a fast-forward of devel.
git push origin devel
git push origin main

# 2. Move the tag. Safe: no Release and no PyPI artifact depend on 3739dfc.
git tag -f -a v1.0.10 -m "LiftOn v1.0.10" <release-sha>
git push --force origin v1.0.10

# 3. Wait for the tests.yml run on that commit to go green. publish.yml
#    deliberately does not re-run the suite, so this is the gate.

# 4. Cut the Release. THIS IS THE IRREVERSIBLE STEP: it fires publish.yml,
#    which uploads to real PyPI via OIDC. A PyPI version can never be
#    re-uploaded, so verify `lifton -V` in the built wheel first.
gh release create v1.0.10 --title "LiftOn v1.0.10" --notes-file <notes>
```

## Follow-ups this release does not cover

- **The live v1.0.10 technical report on khchao.com describes the `3739dfc`
  build.** Its accuracy, completeness and validity claims still hold (see above),
  but the performance panel and the file-size figures moved. It needs a refresh
  and a separate deploy.
- **`docs/source/content/benchmarks.rst`** is still titled "Benchmarks (v1.0.9)".
  Updating it needs benchmark data this release deliberately did not regenerate.
- **The 21.8 % Step-7 SQL collapse** remains reverted. Landing it needs an
  `all_children_full` cache on `MaterialisedLocus`, populated in both walk twins
  and served by `_LFeatureDbProxy`; without it the query signature is uncached
  and every threaded run raises `NotImplementedError` — correct serially, broken
  in parallel, which is exactly what the 24-cell matrix caught.
