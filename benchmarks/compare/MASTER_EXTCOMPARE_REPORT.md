# Broader-comparison sub-study: homology lift-over vs coordinate LiftOver (and CAT)

*Scope: a representative 5-pair, chromosome-subset study added alongside the committed
17-genome 4-way benchmark. The committed `fourway_results.json` is untouched; all new
data lives in `extcompare_results.json` + `/ccb/salz3/kh.chao/lifton_extcompare/`.*

## 1. Why this study

The headline LiftOn v1.0.9 benchmark compares the three **homology lift-over** tools that
take only `(reference genome, reference annotation, target genome)` and emit a re-scorable
target annotation:

- **Liftoff** — DNA-only lift-over (minimap2);
- **LiftOn v1.0.8** — the previous dual-evidence release;
- **LiftOn v1.0.9** — the current release.

miniprot is **not** a peer in that comparison — it is LiftOn's *internal protein evidence*
(the `-M` input LiftOn fuses), so it is shown only as **protein-evidence context**, never as
a headline bar (this demotion already shipped to the technical report).

Two readers asked the natural follow-up: *how does LiftOn compare to tools outside the
"only needs (ref, target)" class?* This sub-study adds two such tools:

- **CAT** (Comparative Annotation Toolkit) v1 and v2 — alignment-based comparative
  annotators that require a **precomputed whole-genome Cactus alignment (HAL)**;
- **UCSC-style coordinate LiftOver** (CrossMap over a chain) — pure coordinate projection,
  **no frame correction**.

Both are **out-of-class**: they need inputs (a Cactus HAL; a chain) that the homology
lift-over tools build internally or do not need at all. The study reports them honestly on
that footing.

## 2. Method

**Representative 5 pairs** were selected, one per divergence tier, at **chromosome-subset** scope
(full-genome Cactus was prohibitive — see §4). The Cactus alignment that LiftOver and CAT
require built only on the close/same-species pairs; the three divergent pairs were **Cactus-
intractable** (§4), so the LiftOver comparison reports the 2 tractable pairs and the
intractability is itself a result.

| pair | tier | ref → target | Cactus alignment |
|---|---|---|---|
| `t1_maize_b73_to_mo17` | same-species | maize B73 → Mo17 | ✅ built (255 MB HAL) |
| `t2_human_to_gorilla` | close | human → gorilla | ✅ built (63 MB HAL) |
| `t3_dog_to_cat` | distant | dog → cat | ❌ cons > 700 GB |
| `arabidopsis_to_rice` | very-distant (plant) | *Arabidopsis* → rice | ❌ cons > 279 GB |
| `t4_human_to_chicken` | very-distant | human → chicken | ❌ (larger/more divergent) |

**Per-pair Cactus alignment.** For each pair we build our own pairwise **Cactus v3.2.1**
alignment (`singularity exec cactus … --consMemory 300G`) from the exact subset FASTAs the
benchmark already extracted. That HAL feeds **both** comparators:
- **LiftOver**: `cactus-hal2chains` → ref→target chain → `CrossMap gff` lift of the reference
  GFF3 → strip CrossMap's leading `chr` from seqids → target GFF3.
- **CAT**: `luigi --module cat RunCat --hal=… --ref-genome=ref --config=…`.

**Scoring.** Every tool is re-scored by the **same tool-neutral parasail kernel** used for the
4-way benchmark (`evaluator.evaluate_tool`), so the numbers are directly comparable. The
driver is `benchmarks/compare/extcompare.py` → `extcompare_results.json`.

**LiftOver — dual metric (the honest out-of-class framing).** Coordinate lift does not
correct frame, so the re-translated CDS scores low. We therefore report LiftOver with **two**
numbers:
1. **re-scored protein identity** — expected LOW (the same parasail score as every tool); and
2. **feature-lift rate** — the fraction of reference coding transcripts whose CDS lifted at
   all (by id). CrossMap fragments a feature at every chain-block boundary, so
   `extcompare._collapse_liftover` first collapses the lifted GFF to one clean transcript per
   reference id (the dominant target seqid, deduped segments) before scoring and before
   counting the lift rate.

miniprot is carried through scoring as the **protein-evidence** series (a diamond marker over
the bars), never a headline bar.

## 3. Results

`extcompare_figures.py` → `figures/extcompare/rfig_extcompare.png` (3 panels: A mean protein
identity, B coding completeness, C the LiftOver dual metric). Numbers below are read directly
from `extcompare_results.json`.

### Mean protein identity (tool-neutral re-score)

| pair | Liftoff | LiftOn v1.0.8 | **LiftOn v1.0.9** | LiftOver | *miniprot (evidence)* |
|---|---|---|---|---|---|
| maize B73→Mo17 (same-sp) | 0.917 | 0.946 | **0.950** | 0.841 | *0.938* |
| human→gorilla (close)    | 0.963 | 0.974 | **0.979** | 0.920 | *0.964* |
| dog→cat (distant)        | — | — | — | — | — |
| *Arabidopsis*→rice (v-distant) | — | — | — | — | — |
| human→chicken (v-distant) | — | — | — | — | — |

*The three distant/very-distant pairs have no LiftOver/CAT column: the per-pair Cactus
alignment both require **could not be built** (cons > 700 GB; see §4). The homology lift-over
columns would be available, but with no LiftOver/CAT to compare against on those pairs they are
omitted from this LiftOver-focused sub-study; they are already covered in the 17-genome 4-way.*

### Coding completeness (recall) and LiftOver feature-lift rate

| pair | Liftoff | v1.0.8 | **v1.0.9** | LiftOver compl. | **LiftOver feature-lift rate** | *miniprot* |
|---|---|---|---|---|---|---|
| maize B73→Mo17 | 0.961 | 0.961 | **0.961** | 0.855 | **0.855** | *0.985* |
| human→gorilla  | 0.999 | 0.999 | **0.999** | 0.994 | **0.994** | *0.999* |
| dog→cat / Arab→rice / human→chicken | — | — | — | — | — | *(Cactus intractable; §4)* |

### What the two confirmed pairs already show

- **LiftOn v1.0.9 leads on protein identity** on both confirmed pairs (0.950 maize, 0.979
  gorilla), above Liftoff (the DNA baseline) and above the miniprot protein evidence, and
  above v1.0.8 — the same ordering the headline benchmark reports, now also true against an
  out-of-class coordinate tool.
- **LiftOver is genuinely out-of-class.** It lifts most features (99.4% of reference coding
  transcripts on the close human→gorilla pair) yet its re-translated CDS scores ~0.06 lower
  than even Liftoff (0.920 vs 0.963) because it does not correct frame. On the
  structurally-variable same-species maize pair its feature-lift rate itself falls to 0.855
  (the Cactus chain simply does not cover the rearranged regions), and identity falls to
  0.841. **This is the honest, informative result: coordinate lift-over answers a different
  question** ("where does this interval land?") than homology annotation lift-over ("what is
  the correct gene model on the target?").
- **The trend the two pairs anchor:** even from same-species → close, LiftOver's quality is
  governed entirely by the underlying Cactus alignment — and that alignment *cannot be built at
  all* on the divergent pairs (§4). So the LiftOver/CAT comparison is structurally available
  only where the genomes are close enough for Cactus to converge in feasible memory — which is
  precisely **not** the regime where annotation lift-over is hardest or most valuable.

## 4. Per-tool feasibility & cost (the honest verdicts)

### CAT (v1 and v2) — **attempted, intractable on this input; no column**

Both CAT versions were **installed and import-clean** (`cat_v2` = CAT v2.2.1 on py3.8 with
Augustus 3.3.3, toil 5.12, the full kent/UCSC toolset; `cat_v1` = CAT v1.0 on py2.7). On the
gorilla HAL, `RunCat` ran **18 of 21 luigi tasks** but its toil sub-workflows
(`Chaining` / `GenerateHints` → `AlignTranscripts` / `Consensus`) **fail opaquely** on the
RefSeq subset — the jobs die without a captured error and toil cleans the jobStore, after a
long sequence of fixes already applied (openssl-3 kent libs via an `LD_LIBRARY_PATH` shim;
an orphan-CDS-free input GFF3; `sqlalchemy<2`; `toil==5.12`; mRNA `gene_biotype`/`gene_id`/
`transcript_id` enrichment). After a substantial time-box this was judged **intractable for
this study** and dropped — `extcompare_results.json` records it as
`cat: { status: "attempted_intractable" }`. This is itself an honest finding: **CAT carries
a heavy operational burden** (a precomputed Cactus alignment *plus* a fragile multi-tool
toil/Augustus/kent stack), in sharp contrast to LiftOn/Liftoff, which need only
`(ref, ref-annotation, target)`.

### Coordinate LiftOver — **working end-to-end**

`cactus-hal2chains` → `CrossMap gff` → seqid-normalize → `_collapse_liftover`. Validated on
human→gorilla (PI 0.920, feature-lift 0.994) and maize (PI 0.841, feature-lift 0.855). The
only non-obvious engineering: CrossMap fragments features and prepends `chr` to seqids; both
are handled in the harness.

### Cactus — **intractable on the divergent subsets (the decisive finding)**

The per-pair Cactus alignment is the bottleneck — and on the divergent pairs it is not merely
slow, it is **infeasible**:

- **Close / same-species pairs build fine.** human→gorilla (63 MB HAL) and maize B73→Mo17
  (255 MB HAL) converged with `cactus_consolidated` inside the default ~90 GiB.
- **Divergent pairs blow up memory super-linearly.** toil auto-estimates the cons-job memory
  from input size (~89.8 GiB here) and **enforces it as an RLIMIT_AS cap**; the divergent
  melting rounds blow past it (`Malloc failed … Cannot allocate memory`) *despite >350 GB
  free*. We confirmed this is a real cactus-graph explosion, not a config slip:
  - `--restart` does **not** help — toil bakes the 89.8 GiB into the saved cons job and
    re-OOMs (verified);
  - a **clean** dog→cat run at `--consMemory 300G` issued the cons job at **279.4 GiB** (so the
    flag applied) and **still OOM'd**, failing even a 24 KB allocation → `cactus_consolidated`
    genuinely exhausted ~279 GB on a **380 MB** subset (a ~700× blow-up; one melting round alone
    destroyed 566 K alignment blocks);
  - it is **not** a masking artefact — the subset FASTAs are 46–49 % soft-masked and Cactus adds
    ~31 % RED masking;
  - a final **clean dog→cat run at `--consMemory 700 GB` (run alone on the 1 TB box) OOM'd again**
    (4 mallocs failed). dog→cat is the *smallest* of the three divergent subsets, so
    *Arabidopsis*→rice and human→chicken (both larger and more divergent — *Arabidopsis*→rice
    already OOM'd at 279 GiB) are at least as infeasible. We stopped there rather than chase
    >700 GB on a shared machine.

**Net:** the Cactus alignment that **both** CAT and chain-based LiftOver depend on is
**operationally intractable for divergent genomes at this scope** — which is the same root
reason CAT could not be evaluated and the LiftOver columns end at the close tier. This is the
sub-study's sharpest result: the entire **alignment-based** comparison class is gated by a
whole-genome alignment that does not scale to divergence, exactly the regime where homology
annotation lift-over is hardest. LiftOn and Liftoff need only `(ref, ref-annotation, target)`
and ran on all 17 whole genomes.

## 5. Take-aways

1. **Within the homology lift-over class, LiftOn v1.0.9 leads** — confirmed again here against a
   coordinate tool, on top of the 17-genome headline result (gorilla 0.979, maize 0.950).
2. **Coordinate LiftOver is complementary, not competitive** for annotation: high feature-lift
   rate (0.994 close / 0.855 same-species), low protein identity (0.920 / 0.841 — no frame
   correction). Reported with both numbers it is an honest reference point, not a rival.
3. **The alignment-based class does not scale to divergence.** Both CAT and chain-based LiftOver
   depend on a per-pair Cactus alignment, and that alignment is **intractable on the divergent
   subsets** (`cactus_consolidated` > 700 GB for a 380 MB dog→cat subset; *Arabidopsis*→rice and
   human→chicken at least as bad). CAT additionally fails in its own toil sub-workflows even on
   the close pair. So the broader comparison is *structurally* confined to the close/same-species
   tier.
4. **The input-class distinction is the real story** — and it is *decisive*, not cosmetic:
   LiftOn/Liftoff need only `(ref, ref-annotation, target)` and ran on all 17 whole genomes
   including the very-distant ones; CAT needs a precomputed whole-genome alignment that
   **cannot be built** for those same divergent pairs; LiftOver needs a chain (from that same
   alignment) and is coordinate-only. The honest comparison is therefore: LiftOn leads where all
   tools can run, and is the *only* class that runs at all where divergence is high.

---
*Reproduce: `python -m benchmarks.compare.extcompare [pairs]` (scoring) ·
`python benchmarks/compare/extcompare_figures.py` (figures). Cactus/LiftOver build scripts:
`/ccb/salz3/kh.chao/lifton_extcompare/bin/{recover_one.sh,make_liftover.sh}`. Nothing here is
committed/deployed without explicit approval.*
