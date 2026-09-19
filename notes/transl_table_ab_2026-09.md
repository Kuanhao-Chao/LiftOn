# `transl_table`: what changed, on whole genomes

Two whole-genome arms, each pinned to its own detached worktree, sharing one
cached `-L`/`-M` so the only difference is the code. The harness asserts which
build each arm loaded, after a first attempt silently ran both arms out of the
main repo — `python -c` puts the working directory at `sys.path[0]`, ahead of
`PYTHONPATH`, so a launch directory containing a `lifton/` package wins and the
comparison proves nothing.

## rice — must not move

The rice reference declares `transl_table=11` on 124 CDS. Table 11 differs from
table 1 only in its start codons, so with `cds=False` it translates
identically and its stop set is the same. The honest prediction is therefore
"no change", and that is a real test of the plumbing rather than a formality:
the run must detect the table, route 96 plastid proteins through it, and still
produce the same bytes.

| | |
|---|---|
| sidecar written | yes — 96 proteins, all table 11 |
| output | **byte-identical**, md5 `1d0e084b2df3f397431d9e626b4746f1` |

So the claim that this fix does nothing to plant plastid genes is measured, not
asserted.

## bee — must move, and only where it should

The bee reference declares `transl_table=5` (invertebrate mitochondrial, where
TGA is tryptophan) on 13 CDS. The target assembly **has no mitochondrion**, so
those genes lift onto nuclear NUMTs — degraded mitochondrial segments in the
nuclear genome. That makes it a harder case than a clean organellar lift, and
the fix still reaches it, because the error was in reading the *reference*.

| transcript | protein identity | mutation |
|---|---|---|
| KEF36_p02 | 0.083 → **0.911** | `stop_codon_gain` gone |
| KEF36_p10 | 0.170 → **0.868** | → `nonsynonymous` |
| KEF36_p01 | 0.212 → **0.752** | → `stop_missing` |
| KEF36_p06 | 0.065 → **0.485** | → `stop_missing` |
| KEF36_p13 | 0.266 → 0.480 | |
| KEF36_p11 | 0.165 → 0.335 | |
| KEF36_p12 | 0.029 → 0.219 | |
| KEF36_p07 | 0.915 → 0.958 | |
| KEF36_p06_1, p06_2, p05_1 | 0.058 → 0.045, 0.100 → 0.092 | down |

Ten improved, three slightly worse. The three are NUMT copies that are
genuinely degraded — their alignments are poor under either code (0.04–0.10),
and scoring them against the right table moves them by a hundredth. Reporting
them is the point: a fix that only ever improves a number is usually a fix that
is measuring itself.

**Containment**: 23,437 mRNA in both arms — 0 added, 0 lost — and **0
transcripts changed outside the 14 that declare table 5**.

## Scope, stated plainly

Affected: annotations declaring a genetic code that reassigns a codon —
mitochondrial tables (2, 5, and the rest), and the alternative nuclear codes.
In this corpus that is human and mouse (13 CDS each, table 2) and bee (13,
table 5).

Not affected: table 11, which is where the plant plastid genes are — 124 in
rice, 104 in arabidopsis. They translate identically to table 1 and are pinned
that way by a test, so the scope of this fix is not overstated.
