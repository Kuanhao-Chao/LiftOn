# Windowed aligner: index the reference over the query's k-mers, not all of them

*2026-09-21. Cycle 3, item P3.*

## Where the time is

Profiling a dog → cat lift put `windowed_align._align_region` at 357 s of a
3,299 s Step 7 (10.8 %), with `_unique_anchors` the largest single component at
116 s **self** time.

Scoped honestly: `windowed_align` is **absent from rice's profile entirely**.
The windowed path only engages above 2,500 aa / 8,000 nt, so this helps
mammalian and divergent lifts and does nothing for a same-species run.

## What it was doing

An anchor is a k-mer that occurs exactly once in the query **and** exactly once
in the reference. The function built a full count-and-first-index map of every
k-mer in each sequence, then intersected.

A k-mer absent from the query — or repeated in it — can never become an anchor.
Indexing it was pure cost, and on a long reference region that is most of the
work.

## The change

Build the query map first and delete its repeats; then scan the reference,
recording only k-mers that survived on the query side.

This is the same anchor set, not an approximation of it. An anchor requires
uniqueness in both sequences, so restricting either side to the other's unique
k-mers leaves the intersection untouched. The sort is unchanged too: an entry's
`q_pos` determines its k-mer (`query[i:i+k]`), so no two anchors share one and
ordering by the tuple is ordering by `q_pos` either way.

## Measured

Alternated arms, anchors asserted equal on every case:

| | length | k | divergence | speedup |
|---|---:|---:|---:|---:|
| protein | 3,000 | 8 | 0.02 / 0.15 | 1.38× / 1.64× |
| protein | 8,000 | 8 | 0.02 / 0.15 | 1.32× / 1.62× |
| protein | 20,000 | 8 | 0.02 / 0.15 | 1.31× / 1.56× |
| DNA | 10,000 | 15 | 0.02 / 0.15 | 1.31× / 1.60× |
| DNA | 50,000 | 15 | 0.02 / 0.15 | 1.53× / 1.82× |
| DNA | 120,000 | 15 | 0.02 / 0.15 | 1.45× / 1.90× |

The margin widens with divergence, because fewer shared unique k-mers means
more of the reference index was being built for nothing — and divergence is the
regime the windowed aligner exists for.

## Exactness

`tests/test_windowed_align.py::TestUniqueAnchorsReferenceRestriction` checks the
implementation against a naive transcription of the original definition, over
inputs chosen to hit each branch (identical sequences, a k-mer repeated on
either side, no shared k-mers, sequences shorter than k, the DNA seed length,
and randomised inputs at four divergences).

The load-bearing gate is a whole-genome dog → cat `cmp`, not a fixture.
