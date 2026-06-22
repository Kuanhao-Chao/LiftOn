# LiftOn Benchmark Pipeline — Design & Implementation

A reproducible, automated pipeline that benchmarks the LiftOn genome-annotation
lift-over tool against the two tools it builds on — **Liftoff** (DNA-level) and
**miniprot** (protein-level) — across six reference→target species pairs, and
quantifies how much LiftOn's protein-maximization fusion improves annotation
quality. This document describes the methodology and implementation; the
companion **`benchmark_comparison_report.{md,html}`** holds the measured results.

---

## 1. Goal & approach

LiftOn lifts an annotation from a reference genome to a target genome by fusing a
DNA-based lift (vendored Liftoff → minimap2) with a protein-based lift (miniprot)
via an in-house CDS-chaining + ORF-rescue algorithm. To measure the value LiftOn
adds, this pipeline runs **all three tools independently** and evaluates each
output with a **uniform metric** so the comparison is apples-to-apples:

- **Liftoff** standalone → `liftoff.gff3`
- **miniprot** standalone → `miniprot.gff3`
- **LiftOn** reusing both via `-L`/`-M` → `lifton.gff3`
  (so LiftOn is exactly the fusion of the same Liftoff + miniprot outputs being
  compared, isolating the contribution of its chaining algorithm).

Each output is scored for three quality dimensions:

1. **Annotation completeness** — fraction of reference (coding) transcripts the
   tool recovered in the target.
2. **Protein-sequence identity** — translated-CDS identity of each lifted
   protein-coding transcript vs the reference protein.
3. **Gene DNA-sequence identity** — transcript-DNA identity of each lifted
   transcript vs the reference transcript.

To keep iteration fast, each benchmark is reduced to **one chromosome pair**.

---

## 2. The six benchmarks

`GRCh38 → CHM13 (MANE)` is the headline benchmark; the reference annotation is
human = **MANE** (canonical one-transcript-per-gene set), all others = RefSeq.

| id | reference → target | type | reference annotation |
|---|---|---|---|
| **`human_mane`** (MAIN) | Human GRCh38 → T2T-CHM13v2.0 | same-species | MANE v1.3 |
| `bee` | *Apis mellifera* HAv3.1 → ASM1932182v1 | same-species | RefSeq |
| `mouse` | *Mus musculus* GRCm39 → NOD_SCID | same-species | RefSeq |
| `human_to_chimp` | Human GRCh38 → chimp NHGRI_mPanTro3-v1.1 | cross-species | RefSeq |
| `mouse_to_rat` | Mouse GRCm39 → rat mRatBN7.2 | cross-species | RefSeq |
| `rice` | *Oryza sativa* IRGSP → ASM3414082v1 | same-species | RefSeq |

Data hub: `/ccb/salz2/jheinz3/shared/lifton/` (genomes, annotations, proteins),
plus CHM13 at `/ccb/salz2/kh.chao/LiftOn_chm13/`, MANE at
`/ccb/salz2/kh.chao/OpenSpliceAI_exp/`, and a ready chr22 subset pair at
`/ccb/salz2/kh.chao/LiftOn/test/`. Exact per-benchmark paths are in
`benchmarks/compare/benchmarks.json`.

---

## 3. Methodology

### 3.1 Subsetting — both genomes to one chromosome (`subset_builder.py`)

Each benchmark is reduced to a single chromosome pair, which shrinks every step
(Liftoff aligns far fewer gene sequences, miniprot far fewer proteins, LiftOn far
fewer loci) while keeping the lift biologically correct.

- **Reference chromosome.** `human_mane` is fixed to `chr22`. For the rest,
  `AUTO_LARGEST_CODING` streams the reference GFF and picks the most mRNA-dense
  real chromosome (skipping unplaced scaffolds).
- **Target chromosome(s).** `AUTO_SYNTENIC` aligns the chosen reference
  chromosome to the **full** target genome with
  `minimap2 -x asm20 --secondary=no`, sums residue-matches per target sequence,
  and keeps the dominant syntenic target chromosome(s) (≥ `max(10 kb, 10 %` of
  total matches`)`). This works for both same-species (1:1) and cross-species
  (where chromosome numbers differ); non-syntenic reference genes legitimately
  count as "lost".
- **FASTA subsetting** uses `samtools faidx`; **GFF subsetting** is a streaming
  filter on column 1 (preserving `##gff-version` / matching `##sequence-region`).
- **Seqid harmonization.** The reference genome header is made to match the GFF
  seqid (asserted: `faidx_seqids ⊇ gff_seqids`). MANE already uses `chr22`, and
  the ready chr22 FASTAs are already `>chr22`, so the MAIN benchmark needs no
  rewrite.
- **Reference proteins for the subset.** For `human_mane` (`transcript` target
  space) the CDS are translated from the reference genome with headers = the
  mRNA id (`rna-NM_…`) so miniprot's `Target` maps 1:1. For the others
  (`protein` space) the provided protein FASTA (NP_/XP_ accessions) is filtered
  to the subset's protein accessions.

Outputs land under `work/<id>/subset/` with a `subset.manifest.json` recording
the chosen chromosomes, feature counts, header maps, and the
`protein_accession → mRNA_id` table.

### 3.2 Running the three tools (`tool_runners.py`)

Every invocation is wrapped by `run_profiled` (reused from the Phase-16 harness,
`benchmarks/run_benchmarks.py`), which runs it under `/usr/bin/time -v` and
records wall-clock + peak RSS. Each subprocess gets a composed `PATH`
(`/ccb/sw/bin:/home/kh.chao/bin` prepended) and `PYTHONNOUSERSITE=1`. Stages are
guarded by `.done/<tool>.done` sentinels for resumability.

```bash
# Liftoff (standalone)
liftoff -g ref.chrom.gff3 -o liftoff.gff3 -u unmapped.txt -dir intermediate \
        -p 8 -copies -polish  tgt.chrom.fa  ref.chrom.fa
# (-polish appends _polished to the output name; the runner normalizes it.)

# miniprot (standalone)  — --gff-only writes GFF3 to stdout
miniprot -t 8 --gff-only  tgt.chrom.fa  ref.proteins.subset.faa  > miniprot.gff3

# LiftOn (reuses both via -L / -M)
lifton -t 8 -copies -ad RefSeq -g ref.chrom.gff3 \
       -L liftoff.gff3 -M miniprot.gff3 -o lifton.gff3 \
       tgt.chrom.fa  ref.chrom.fa
```

LiftOn is **not** run with `-E` (its built-in evaluation mode has a known
`mapped:0` parsing bug); all metrics come from the custom evaluator below.

### 3.3 Evaluation — custom parasail-kernel evaluator (`evaluator.py`)

The primary metric reuses **LiftOn's own kernel** so the identity numbers match
LiftOn's internal definitions exactly:

- **Sequence extraction** uses `lifton.extract_sequence.get_dna_sequence` /
  `get_protein_sequence` (strand-aware exon/CDS splicing, N-padding, translation)
  — applied identically to the reference (from the reference genome) and to each
  tool's lifted features (from the **target** genome).
- **Identity** uses `lifton.align.trans_align` / `protein_align` (parasail
  Needleman–Wunsch: DNA `ACGTN*` matrix gap 5/2; protein BLOSUM62 gap 11/1) +
  `lifton.get_id_fraction` (gap-compressed BLAST identity).
- **Per-tool ID mapping** (`id_mapping.py`): Liftoff/LiftOn carry the reference
  transcript id (a trailing `_N` copy suffix is stripped only if the base exists
  in the reference set); miniprot's `Target` attribute carries the reference
  protein accession (mapped to the mRNA via the subset GFF) or, for MANE, the
  mRNA id directly. Multiple lifted copies of one reference transcript keep the
  max-identity copy as "recovered"; the rest count as `extra_copies` and never
  inflate completeness.

**Metric definitions.**

- `completeness_coding` = recovered coding transcripts / reference coding
  transcripts. (`completeness_all` additionally counts non-coding for
  Liftoff/LiftOn.)
- Protein identity is CDS-translated and **directly comparable across all three
  tools**.
- DNA identity basis differs by tool and is flagged in every output: miniprot has
  no exon/UTR records, so its DNA identity is **CDS-spliced** (`dna_basis="cds"`);
  Liftoff/LiftOn use the **full transcript** (exon-spliced, `dna_basis="transcript"`).
  Reference DNA uses the matching basis per tool so each comparison is internally
  consistent.

Outputs per tool: `<tool>.transcripts.tsv` (per-reference-transcript row) and
`<tool>.summary.json` (completeness, identity mean/median/%-identical + 10-bin
histogram, extra copies, joined wall-clock/peak-RSS).

### 3.4 Cross-check — LiftoffTools `variants` (`liftofftools_wrap.py`)

As an independent confirmation, `liftofftools variants` is run on the Liftoff and
LiftOn outputs (miniprot is excluded — its GFF lacks the gene→mRNA→CDS hierarchy
`variants` requires). Its `variant_effects` classes (`identical`, `synonymous`,
`nonsynonymous`, `frameshift`, `start_lost`, …) are parsed and reported beside
the custom metric. Note the two are correlated but not identical: LiftoffTools
`identical` means *no variant at all* (DNA **and** protein identical), whereas the
custom `%identical` counts protein identity = 1.0 (which includes synonymous DNA
changes), so the custom count is expected to be ≥ LiftoffTools `identical`.

---

## 4. Implementation architecture

All pipeline code lives under `benchmarks/compare/` (a sub-package that reuses the
LiftOn library and the Phase-16 profiling helpers):

| Module | Responsibility |
|---|---|
| `run_compare.py` | Driver/CLI; per-benchmark stage orchestration; report trigger |
| `benchmarks.json` | Registry: tool binaries + per-benchmark data paths & chrom policy |
| `subset_builder.py` | Chromosome subsetting, `AUTO_SYNTENIC` minimap2 detection, seqid harmonization, subset proteins |
| `tool_runners.py` | Profiled standalone Liftoff/miniprot + LiftOn-reusing-both runners |
| `id_mapping.py` | Per-tool lifted-mRNA → reference-mRNA resolvers |
| `evaluator.py` | Custom completeness + protein/DNA identity (reuses LiftOn kernel) |
| `liftofftools_wrap.py` | LiftoffTools `variants` cross-check + parser |
| `reporter.py` | Comparison tables + matplotlib plots (base64-inlined) + markdown; HTML via `build_spec_html.py` |
| `profiling.py` | Re-exports the audited `/usr/bin/time` parsing from `run_benchmarks.py` |

Working data is written under `benchmarks/compare/work/<id>/{subset,tools,eval}/`,
and the final deliverables to `plans/` (`benchmark_comparison_report.{md,html}`,
`benchmark_comparison.tsv`, `benchmark_plots/`).

---

## 5. How to run

```bash
conda activate lifton_devel        # the LiftOn dev env (has lifton + a working parasail)
cd /ccb/salz3/kh.chao/LiftOn
export PYTHONNOUSERSITE=1

# Everything (resumable; ~minutes per benchmark on the chromosome subsets):
python -m benchmarks.compare.run_compare --all -t 8

# A single benchmark, or a single stage:
python -m benchmarks.compare.run_compare --benchmarks human_mane
python -m benchmarks.compare.run_compare --benchmarks bee --only-subset
python -m benchmarks.compare.run_compare --all --only-eval     # re-score + rebuild report
```

Every stage writes a `.done/<stage>.done` sentinel, so re-runs skip completed
work; `--force` ignores the sentinels.

### Environment note

The system/base-env `liftoff` and the PATH `liftofftools` are broken on a
numpy/parasail incompatibility (`numpy has no attribute 'intc'`). The pipeline
therefore uses a standalone Liftoff installed into `lifton_devel`
(`pip install --no-deps Liftoff`, which reuses that env's working parasail) and
the working `liftofftools` **env** binary
(`/home/kh.chao/miniconda3/envs/liftofftools/bin/liftofftools`). Both paths are
pinned in `benchmarks.json`.

---

## 6. Results

The measured comparison (completeness, protein identity, DNA identity, runtime,
peak RSS for Liftoff vs miniprot vs LiftOn across all benchmarks) is generated
into **`plans/benchmark_comparison_report.html`** (self-contained, with plots)
and **`plans/benchmark_comparison.tsv`**. All six benchmarks completed; mean
**protein-sequence identity** (the headline metric) per benchmark:

| Benchmark | Type | Completeness (coding) | Liftoff | miniprot | **LiftOn** | LiftOn DNA id | Δ protein (LiftOn−Liftoff) |
|---|---|---|---|---|---|---|---|
| **`human_mane`** (MAIN) | same-sp | 99.8 % (422/423) | 0.9931 | 0.9897 | **0.9942** | 0.9979 | **+0.001** |
| `bee` | same-sp | 99.7 % (3122/3130) | 0.9825 | 0.9807 | **0.9916** | 0.9964 | **+0.009** |
| `mouse` | same-sp | 94.0 % (7677/8170) | 0.9903 | 0.9796 | **0.9941** | 0.9939 | **+0.004** |
| `human_to_chimp` | cross-sp | 99.0 % (12718/12842) | 0.9649 | 0.9654 | **0.9806** | 0.9842 | **+0.016** |
| `mouse_to_rat` | cross-sp | 90.4 % (7384/8170) | 0.8588 | 0.8763 | **0.9104** | 0.8748 | **+0.052** |
| `rice` | same-sp | 100 % (5850/5850) | 0.9986 | 0.9952 | **0.9990** | 0.9999 | **+0.000** |

(Protein identities are means over recovered coding transcripts; miniprot DNA
identity is CDS-based. Completeness/DNA identity are shared between Liftoff and
LiftOn because LiftOn reuses Liftoff's exon/transcript structure via `-L` — its
contribution is the protein-sequence refinement.)

**LiftOn's mean protein identity meets or exceeds Liftoff's (and miniprot's) on
every benchmark**, and the gain grows with evolutionary distance: negligible on
near-identical same-species pairs (rice, human MANE), but **+0.016 human→chimp
and +0.052 mouse→rat**. The protein-maximization fusion recovers protein-sequence
accuracy that the DNA-only lift misses — most valuably in the cross-species regime
where Liftoff's nucleotide alignment degrades — without sacrificing completeness.
Full per-tool tables, identity distributions, resource usage, and the LiftoffTools
`variants` cross-check are in the comparison report.

---

## 7. Design choices & limitations

- **Chromosome subset** trades genome-wide coverage for speed (per the project
  goal of fast iteration). Cross-species benchmarks subset both genomes to one
  syntenic chromosome, so reference genes whose ortholog lies on a different
  target chromosome count as lost — expected, and reported via the loss columns.
- **miniprot is protein/CDS-only** (no UTR, no non-coding): its DNA identity is
  CDS-based and its completeness denominator is coding-only; non-coding reference
  transcripts are unrecoverable by miniprot by construction. Flagged throughout.
- **Apples-to-apples LiftOn**: by reusing the exact standalone Liftoff + miniprot
  outputs (`-L`/`-M`), the comparison isolates LiftOn's chaining contribution
  rather than conflating it with version differences in its internal aligners.
- The evaluator reuses LiftOn's **own** alignment/identity kernel, so the
  reported identities are defined identically to LiftOn's internal metrics; the
  independent LiftoffTools cross-check guards against a metric-definition bias.
