# Open-issue triage, re-verified against `devel` for the v1.0.12 release

Every claim below was checked against the current tree, not against the earlier
triage pass. Two corrections to that pass: PyPI's latest is **v1.0.11**
(v1.0.10 was tagged but never published), and anything described as v1.0.12 is
unreleased until the tag lands — so these replies are to be posted **after** the
release, not before.

| # | Verdict | Action |
|---|---|---|
| 26 | already fixed | close |
| 25 | already fixed | close |
| 37 | fixed in v1.0.12 | reply, close on confirmation |
| 14 | half already fixed, half fixed in v1.0.12 | reply, close on confirmation |
| 23 | needs the reporter | reply, keep open |
| 16 | tool fixed, published file not | keep open (ops) |
| 38 | half fixed, half open | keep open (enhancement) |
| 31 | not implemented | keep open (enhancement) |

---

## #26 — Duplicate exons, overlapping CDSs in output — CLOSE

The duplicate exon rows and doubled CDS came from the Case-1 branch of
`update_cds_list`: the guard that stops an already-processed, CDS-bearing exon
from being re-emitted was gated on `optimize`, so with `optimize=False` every
downstream exon emitted another copy of it. Fixed in `7e3093b`, with
`tests/test_update_cds_list_case1.py` covering both `optimize` settings by
asserting the emitted exon list rather than only the identity score.

Three later changes make this class of output impossible to ship at all: an exon
is extended to cover its CDS and the exons are sorted at the write funnel; exon
IDs are renumbered 5′→3′ whenever they would collide; and structural validation
now runs before publication, where a duplicate ID or an out-of-order coordinate
is an error that blocks the write. So a run that hit this today would fail
loudly rather than emit a file gffread cannot read.

Please upgrade to v1.0.11 or later. Reopen if you still see it.

## #25 — `No lines parsed` and `'LiftOn_FEATURE' object has no attribute 'ref_gene_id'` — CLOSE

Both halves are fixed.

The `No lines parsed` pair came from a reference sequence-name mismatch that
produced empty FASTAs. LiftOn now warns per missing sequence name while
extracting, and raises a named `LiftOnInputError` when features were iterated
but both FASTAs came out empty — before either aligner runs. A miniprot that
errors or returns nothing no longer produces two cryptic database failures
either; the run degrades to a Liftoff-only result.

The `AttributeError` had a specific cause worth recording: a three-level
reference hierarchy — `gene → primary_transcript → miRNA → exon`, common in
RefSeq — put a generic feature where a gene was expected. That hierarchy is now
handled directly, and the write gate reads the attribute defensively, so it
cannot abort a run. This is why it hit 5 of your 12 genomes and not the rest.

Please upgrade to v1.0.11 or later.

## #37 — `GFF does not contain any gene features` on a bacterial annotation — FIXED

You were right that this was not a problem with your file. Two things were
wrong on our side.

LiftOn auto-detects which top-level types to lift by looking for a type that has
a top-level instance with children. A flat annotation — bakta output, a miniprot
GFF — has top-level `CDS` rows and no `gene` at all, so nothing qualified, and
the detection fell back to `["gene"]`, which selected nothing. It now falls back
to the top-level types the annotation actually has, skipping the ones that
describe a sequence rather than a feature on it (`region`, `chromosome`, …), so
a flat annotation lifts its `CDS` rows.

Second, when a selection really is empty the run used to continue for several
steps and then die inside the vendored Liftoff with the bare message you saw.
It now stops immediately with an error that names the file, what was looked for,
and the feature types the annotation actually contains.

Both land in v1.0.12.

## #14 — Race conditions in parallel runs — MOSTLY FIXED

Better news than the earlier reply gave. Of the two races you reported, the
gffutils one is already fixed, and a third you did not name is too.

- **The shared database.** The reference database is now content-addressed with
  a manifest and a lock file. A run whose manifest matches opens the existing
  database read-only and never rebuilds it, and a build goes to a temporary file
  published with an atomic rename. Pre-building the databases once and then
  fanning out works exactly as you asked.
- **The minimap2 index** (not in your report, but it bites the same workflow).
  The `.mmi` is no longer written next to the input FASTA; it is built inside
  the run's own directory and installed atomically. An index next to the input
  is only ever read.
- **The output file.** The final GFF3 is staged and published atomically, with
  the staging name keyed to the destination.

The race you actually reported — every job with a shared output parent writing
into one `lifton_output/` — was still there. v1.0.12 adds
`-dir/--intermediate-dir`, which gives a run its own directory for intermediate
files, statistics, the score table and the run manifest. Your
`output/$SOURCE/$TARGET.gff` layout becomes:

```
lifton -g ref.gff3 ref.fa target.fa \
       -o output/$SOURCE/$TARGET.gff \
       --intermediate-dir scratch/$SOURCE.$TARGET
```

The default is unchanged, so nothing existing moves.

## #23 — Reference database construction gets stuck — NEEDS A RE-TEST

We could not reproduce a hang, and the evidence points at storage rather than a
loop: the leading process was in uninterruptible sleep, which is where a process
waits on I/O.

Several things have changed that bear directly on what you saw. The build now
prints which file it is building and with which settings before it starts, so a
slow build no longer looks like a silent hang. It no longer blind-retries: only
a recognised duplicate-identifier failure falls through to a second strategy,
where before any failure re-ran the whole multi-minute build up to three times.
The result is cached and published atomically, so an interrupted build cannot
leave a half-written database that stalls the next run.

Could you re-test on v1.0.11 with the database on local disk rather than network
storage, and report the wall time? If it is still slow we will take it from
there.

## #16 — Published CHM13 file has incorrect exons — OPEN (ours to do)

Confirmed, and the tool bug behind it is fixed: colliding exon IDs within a
transcript are renumbered 5′→3′, and structural validation now blocks
publication on a duplicate ID, so a fresh lift cannot produce this file shape.

The published `JHU_LiftOn_v1.0_chm13v2.0.gff3` is unchanged, though — closing
this needs us to re-run the lift on the current release and re-upload. Tracking
it here until that is done. Thank you for the precise report; NOC2L was exactly
the right example.

## #38 — UTRs and mature peptides are not all lifted — PARTLY FIXED

The half about your viral annotation is fixed. A gene with no exons of its own —
`gene → CDS-polyprotein → mature_protein_region_of_CDS` — now keeps its whole
subtree, where previously the subtree was lost. The mapping error you saw from
miniprot also no longer aborts the run; it degrades to a Liftoff-only result.
Worth a re-test on v1.0.11.

The other half is a genuine gap and this issue stays open for it: for an
ordinary `mRNA → exon/CDS/five_prime_UTR` transcript, the UTR rows are dropped.
LiftOn rebuilds a transcript's structure from its exons and CDS when it merges
the two alignments or rescues an ORF, and a UTR's coordinates do not survive
that rebuild, so carrying them through needs the rebuild to recompute them
rather than a switch we can flip. Separately, a UTR row with no parent at all is
not in the lift set, which is a smaller fix.

## #31 — Wrong paralog wins a contested locus — OPEN (enhancement)

Your diagnosis is right, and gene order is the right idea. LiftOn has no synteny
or collinearity model today: each locus is decided on alignment evidence alone,
so where a short 100 %-identity paralog and the correct longer gene compete,
the paralog can win and the correct gene is reported unmapped.

Neither rescue pass covers your case. The miniprot-only rescue fills only loci
no gene occupies, and your candidate overlaps the emitted model by about 19 %,
above the threshold. The cross-locus pass requires the better hit to be on a
*different* sequence; both of yours are on `CITME_001`.

We measured how large this class is while preparing v1.0.12: on distant
transfers it is now 79-91 % of everything still missed, which makes it the main
remaining accuracy problem rather than a corner case. That does not make it
quick — it needs an anchor or gene-order model feeding locus arbitration — but
it is the thing worth building next, and this issue is the reference for it.
