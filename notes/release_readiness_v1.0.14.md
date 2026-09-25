# Release readiness — v1.0.14

*Rewritten 2026-09-23 from the final qualification; updated 2026-09-24 after a
second pre-release scan, a re-check on the reporter's own command, and a
performance campaign against v1.0.13. Supersedes the cycle-3 version of this
note, which is kept as Appendix A (verbatim; headings demoted one level).*

**Verdict: ready to release** — every requalification gate below passed on
the frozen commit; CI on the pushed head: pending (runs on push). Frozen commit **`3ab51cf`** on `v1014-integration` (it
was `e688ff5`; the section below says what changed and how it was
requalified). Nothing has been merged to `main`, tagged, released, uploaded to
PyPI, changed on Bioconda, or published for CHM13, and no reply has been sent.
The runbook at the end lists every one of those steps; each is an outward
action that needs sign-off.

v1.0.14 is 77 commits on top of v1.0.13 (`b2fe59f`, which is `main` and
`devel`), a clean fast-forward. It was built across four Claude Code cycles
(2026-09-14 → 09-21) and one Codex session (2026-09-22/23, "Improve LiftOn
comprehensively"), then reconciled, reviewed, scanned twice and qualified here.

## Second pre-release scan, requalification and performance (2026-09-24)

### What the scan changed

Full findings, with the reproduction for each:
`notes/lifton_v1.0.14_scan_2026-09-24.md`. Every fix has a test that fails
before it and passes after; each is inert where its defect is absent (24-cell
matrix green, no golden edit).

- **A crash on stock NCBI input** (C1): the new reference normalization
  aborted the whole run on a sparse coding model it could not rebuild
  unambiguously — the repo's own GenBank yeast R64 annotation (47 Ty
  gag-pol frameshifts) exited 2 with no output. Such a model is now lifted as
  written (v1.0.13's treatment), warned once and counted; `--strict-gff`
  keeps it fatal.
- **An unsplittable CDS cost its whole gene** in Step 7 and aborted `-E`
  (H1); it now costs the transcript, still counted.
- **Selenocysteine for GENCODE/Ensembl users** (M2): the `transl_except` fix
  read only RefSeq's qualifier. `Selenocysteine` GTF rows (which gffread
  drops) and `stop_codon_redefined_as_selenocysteine` GFF3 rows are now read
  the same way. GENCODE v49 → CHM13: 71 selenoprotein transcripts 0.66–0.68 →
  0.995–0.997, every written codon a TGA, every non-Sec gene unchanged.
- **`transl_except` robustness** (M3, M4, L1, L3, L10, L11): the ORF search
  reads through declared codons; one malformed value no longer discards a
  transcript's others; a value that is neither a read-through nor a start is
  written only onto an identical codon; absurd ranges are refused.
- **Declared genetic code for miniprot-derived models** (M1): mitochondrial
  genes placed from miniprot hits translated with the standard code. dog →
  cat COX2 0.259 → 0.965, CYTB 0.063 → 0.889; human → zebrafish COX1 0.053 →
  0.848.
- **Accounting and speed** (L4, L5, M5): drop counts no longer depend on
  `--threads`; the second-locus sub-pass fell from 24–46 s to under 1 s per
  genome at `-t 1`.
- **Validator** (L7, N2): overlap checked within one sequence and strand, a
  declared ribosomal-slippage overlap is a warning, and the phase check reads
  a 5′-partial CDS's own first phase — without that, v1.0.14's correct output
  drew 5,648 spurious warnings on dog → cat.

On the qualification corpus, `e688ff5` → `3ab51cf` changes only transcripts
that declare a `transl_except` or a genetic code — rice 1, bee 4, dog → cat
13, human → zebrafish 9, drosophila 0 (at `-t 1` and `-t 8`), none worse,
nothing added or removed. CHM13 (1,237,341,536 bytes) and MANE v1.5 → CHM13
are byte-identical to the `e688ff5` outputs.

One fix was written and withdrawn: `1403433` kept a transcript whose two
Liftoff CDS rows share a base, the full suite failed it (the rule it broke —
never count coding bases twice — is pinned by a test), and `1081a2f` reverted
it. It was committed after running only the targeted test file.

### Requalification (`3ab51cf`)

`3ab51cf` differs from the fully requalified `ac7e863` only in the validator's
phase check, which cannot change lift output: a subset lifted with each is
byte-identical (13,766,858 bytes).

| gate | result |
|---|---|
| full suite, Python 3.10 | **2,522 passed**, 14 skipped, 0 failed (49:04) — `ac7e863`'s 2,518 plus the 4 new phase tests |
| full suite, Python 3.11 | **2,522 passed**, 14 skipped, 0 failed (48:56) |
| full suite, Python 3.12 (environment rebuilt under `/ccb/salz3`; the `/tmp` one lost files to tmp cleanup mid-qualification) | **2,522 passed**, 14 skipped, 0 failed (46:33) |
| `make test-fast` (24-cell byte-identity + integration) | 30 passed, no golden edit |
| fatal flake8 | 0 |
| Sphinx docs | succeeds; 77 warnings, the identical set to `ac7e863`'s |
| packaging (wheel/sdist, compiler-free installs, installed chr22 lifts) | `twine check --strict` pass; wheel `lifton/` = tree (95 files, 0 differ; 12 not shipped are vendored Liftoff's tests); compiler-free wheel installs on 3.10/3.11/3.12 and sdist on 3.11; 4 of 4 installed chr22 lifts exit 0, `gff3-validate` exit 0, **byte-identical to each other and to `ac7e863`'s and `e688ff5`'s** (md5 `4d98f5c8…`). wheel sha256 `2ea59cc6…110f`, sdist `36e38ea1…6365` |
| `make benchmark-gate` (isolated export, own `work/human_mane`) | **GATE PASS**: human_mane protein identity 0.99425 → 0.99542, completeness 0.99764 unchanged (the same numbers as `e688ff5`) |
| eight-cell second-locus ladder (`ac7e863`, lift output unchanged since) | SAFETY 8/8, identical line for line to `e688ff5`'s |
| CHM13 regeneration (`ac7e863`) | byte-identical to the staged file; `3ab51cf`'s validator: 0 errors, 11,975 warnings (the 56 phase warnings were N2's false positive) |
| `-t 1` = `-t 8` (rice, drosophila, dog → cat) and `-t 1` = `-t 4` (34 subsets) | byte-identical, all |
| CI on the pushed head | pending (runs on push) |

### The reporter's case, re-checked

The user who reported the selenoprotein problem lifted MANE v1.5 to CHM13 with
`lifton -g MANEv1.5.gff -chroms chrom_mapping.txt -copies -sc 0.9
chm13v2.0.fa hg38.p14.fa`. Re-run exactly that way (fresh) and on the shared
aligner output (cached):

- SEPHS2 is the correct model: CDS chr16:30,830,620–30,831,966, mRNA
  30,829,870–30,832,113 with both UTRs, identity 1.000,
  `transl_except=(pos:complement(30831787..30831789),aa:Sec)` — a TGA in
  CHM13.
- All 25 selenoprotein transcripts: 0.654 → 0.999 (18 below 0.9 → 0; lowest
  GPX1 0.990, whose alanine repeat is shorter in CHM13). 9 non-AUG starts keep
  their scores and gain CHM13 coordinates.
- Independent check, no LiftOn code: each lifted CDS spliced from CHM13 and
  translated with the written `transl_except` applied, against NCBI's own
  GenPept records: 27 of 33 identical; 5 differ where CHM13 and GRCh38
  genuinely differ, 1 where NCBI's protein differs from the GRCh38 genome.
  v1.0.13: 0 of 33.
- v1.0.13 → v1.0.14 on that input changes 40 of 19,898 transcripts: the 34
  declaring ones, DDTL 0.193 → 0.296 and TMEM255B 0.829 → 0.997 (merges that
  now succeed), TMEM52 (an overlapping exon removed), SIRPB1 and STOX1 (the
  CDS score column only), and adds one TEX28 copy at 0.990.
- `lifton_output/liftoff/liftoff.gff3` is Liftoff's intermediate and keeps
  GRCh38 coordinates; the draft reply says so
  (`notes/reply_simon_transl_except.md`, not sent).

### Performance against v1.0.13

`notes/v1.0.14_performance_evaluation.md` (8 whole genomes, 34 subsets, every
tool scored by one evaluator, same aligner input). In short:

- mean protein identity up on all 8 genomes and 32 of 34 subsets (none down);
  paired, 3,484 transcripts better and 607 worse on the genomes;
- +1,068 transcripts on the genomes against 15 given up, every one of the 15
  traced (mostly genes nested in the intron of a neighbour whose model is now
  complete); 962 of the gains are second-locus placements;
- `gff3-validate` errors 21–1,693 per v1.0.13 output → 0 on every v1.0.14
  output; v1.0.13 finished `partial_success` on 4 genomes, v1.0.14 `success`
  on all 8;
- 1.055× faster by geometric mean (0.97–1.20×), peak RSS −1.1 to −1.4 GiB on
  the two largest distant lifts.

## What v1.0.14 changes

### Defaults that change the annotation (each has an opt-out)

| change | opt-out | evidence |
|---|---|---|
| **Second-locus rescue** — place a reference gene at a second target locus when miniprot finds it there and no emitted model reaches it (what a whole-genome duplication leaves) | `--no-rescue-second-locus` | human → zebrafish: 689 genes added, +690 target genes covered by zebrafish's own GRCz11 annotation, 0 lost; rice → sorghum: 51 of 53 placements supported by sorghum's own RefSeq, p = 0.001 (below) |
| **Declared genetic code** (`transl_table`) honoured in extraction, translation, ORF search and stop completion | — (table 1 output unchanged) | the published CHM13 file truncated 4 of 13 human mitochondrial genes; now all 13 match reference length (`notes/transl_table_ab_2026-09.md`) |

### `transl_except` (reported by a user; this round)

RefSeq declares every selenocysteine as `transl_except=(pos:…,aa:Sec)`, and
LiftOn ignored the qualifier. Both the reference and the lifted protein carry
`*` at the UGA, and three decision points read it as a premature stop — the
identity score (an identical selenoprotein scored residue/length: SEPHS2
60/449), the chaining chunk score, and `find_variants` (`stop_codon_gain` →
ORF search) — so the correct Liftoff model lost to a truncated one. On
GRCh38 → CHM13 **all 25 human selenoprotein genes** were affected: 53
transcripts at mean identity 0.666, 36 below 0.9; SEPHS2 lost 118 N-terminal
residues and both UTRs. The carried attribute kept GRCh38 coordinates.

v1.0.14 places each declared codon on the reference protein at Step 3, reads
declared recoded stops (Sec, Pyl, readthrough `Other`, an amino acid over a
stop) through in every scorer, the variant call, chaining and stop
completion, and rewrites the attribute into the lifted model's coordinates at
the one output funnel — dropped where the target codon no longer needs it,
never left in reference coordinates. The benchmark evaluator applies the same
rule to every tool. Measured (details under *Qualification evidence*):

- **CHM13:** all 53 selenoprotein transcripts 0.666 → **0.998** (36 below 0.9
  → 0); SEPHS2 is Liftoff's full model at 1.000 with its Sec at CHM13
  coordinates (a TGA); 1,671 rows changed, all in the 99 declaring
  transcripts; 127 of 127 values placed correctly; nothing worse.
- **Five whole-genome pairs** against the previous build: every changed row
  belongs to a transcript whose reference declares an exception (0 outside);
  every emitted value lies in its CDS, in frame, and for Sec/Pyl/Other on a
  stop codon in the target genome (881 of 881); no declared model got worse
  except one low-quality dog → cat model (−0.007, the ORF rescue's 1 %
  threshold — see limits).
- It is not only selenoproteins: RefSeq also declares stop readthrough and
  genome-error stops this way — drosophila's readthrough isoforms improve 482
  of 486 (0.852 → 0.965), rice's 101 of 107 (0.660 → 0.987).
- **Recall:** the miniprot-only rescue scores with the same rule, so
  human → zebrafish now places **seven selenoprotein genes it used to miss**
  (GPX4, DIO1, DIO2, DIO3, SEPHS2, SELENOT, SELENOM), each on its zebrafish
  ortholog in GRCz11.

### Output-corrective fixes (defects that shipped in earlier releases)

| fix | reach on the corpus |
|---|---|
| `--threads 1` (the default) now equals `--threads N` | rice: 17 organellar genes had duplicated exons and a doubled CDS at `-t 1` |
| no overlapping exons or CDS within a transcript (#26's open half) | 27 / 96 / 47 / 30 / 13 transcripts on CHM13 / h2z / drosophila / rice / bee → 0 |
| an intron-spanning CDS is split into exonic segments with per-segment phase (Codex) — supersedes cycle 4's fix, which was wrong | not reached by the corpus; pinned by serialize-and-validate tests |
| trans-spliced fragments sharing an ID bind their transcripts to the fragment that contains them — cross-sequence (Codex) and same-sequence (this round) | rice `nad5`; drosophila `mod(mdg4)` (5 fragments, 31 transcripts) |
| organellar genes with exons listed twice are indexed for the rescue | 12 rice genes whose rescue candidates were silently abandoned |
| a CDS naming an undeclared `Parent` no longer aborts a lift | 111 such rows in `NCBI_RefSeq_no_rRNA.gff` |
| a worker pool that cannot fork falls back to in-process work | strict-overcommit hosts |
| a transcript whose lifted CDS encodes no protein no longer costs its whole gene (this round; **every release since v1.0.9**) | dog → cat: a 37-transcript gene over one 3-bp CDS |
| a gene's tRNA/rRNA children are not revisited as loci of their own (this round) | dog → cat: 395 false pipeline failures → 0; status `partial_success` → `success`; annotation unchanged |
| `LIFTON_RESCUE_ISOFORM_WORKERS` + a rescue with no other isoform no longer aborts the run (this round; **shipped in v1.0.12 and v1.0.13**) | reproduced on the drosophila ladder cell at `-t 1` |
| a miniprot candidate whose CDS cannot be split is skipped, not allowed to drop the Liftoff gene (this round) | never observed; guarded |
| GTF conversion leaves nothing in the system temp dir (this round) | 78 directories / 99 MB had accumulated from tests |

### Validator (`gff3-validate`, `--validate-output`) — user-visible

New ERROR checks: overlapping exons, overlapping CDS (cycle 3), and **a CDS
outside every exon of its transcript** (`cds_exon_containment`, Codex). A file
that validated before can now fail. Measured before shipping (below): the
containment rule fires on **none** of six RefSeq reference annotations (no
false positives) and on none of the v1.0.14 outputs.

### Accounting, performance, robustness

- Every class of dropped reference feature is counted, reported once, and
  recorded in `run_manifest.json`; the miniprot rescue now counts what it
  abandons (0 on this corpus is now evidence, not silence).
- Step 7 roughly a tenth faster (translation 2.11×, attribute encoding 3.21×,
  cloning 1.70×; byte-identical); windowed-aligner anchors 1.3–1.9×;
  dog → cat Step-7 dispatch −13.7 % (`notes/step7_profile_2026-09-21.md`).
- `--rescue-max-inflight` bounds isoform-rescue memory (−24.7 % peak on h2z
  for +4.4 % wall).
- Sparse coding references and GTF input handled more reliably; alternate-locus
  and patch contigs reported.
- Packaging: `mappy` is an optional extra (v1.0.13, issue #78); no dependency
  changed since v1.0.13.

## This round: reconciling two sessions

The Codex session implemented its v1.0.14 plan in an **uncommitted worktree
in `/tmp`** and ended mid-flight. Its work was secured first (patch +
untracked files, byte-verified, and its 917 MB of evidence copied to
`lifton_improve/v1014_evidence/`), then reviewed change by change.

**Kept as written:** the CDS split, the validator rule, the second-locus
truth tool, the docs/version metadata.

**Changed in review:**

- *An ambiguous trans-spliced family aborted the whole run* (`LiftOnInputError`)
  on user-supplied `-L` input — one gene losing a genome, the Iteration-21
  shape. It is now left as bound, warned once and counted
  (`liftoff_same_seqid_parent_ambiguous`).
- The drop-ledger text and the silent-losses note still described cycle 4's
  wrong behaviour; corrected.
- The new benchmark script was unregistered in the inventory, which would have
  turned `test_benchmark_inventory` red on commit.

**Found by qualifying, fixed here** — each with a test that fails without
it (the ladder's harness fix excepted):

- **A gene lost over one unscorable transcript.** dog → cat reported
  `partial_success` with 396 pipeline failures, identical in the cycle-4
  baseline. One was real: an empty lifted protein (a 3-bp CDS at phase 2)
  reached the aligner's guard, whose error no caller catches, and took its
  37-transcript gene with it — in every release since v1.0.9.
- **395 false failures.** The other 395 were tRNA/rRNA children re-enumerated
  as loci because those types have top-level instances — Iteration 20's
  detect-from-top-level shape at the one site that fix missed. All 395 were
  already in the output; they hid the real loss. Reach, measured before
  either fix: the other five genomes' manifests record no failures at all,
  and the re-enumeration predicate selects exactly these 395 on dog → cat and
  nothing on the other five Liftoff inputs.
- **The isoform-rescue empty-batch crash** — the ladder's first cell died with
  `range() arg 3 must not be zero`. Released since v1.0.12, reachable whenever
  `LIFTON_RESCUE_ISOFORM_WORKERS` is set, which LiftOn's own fork-failure
  warning recommends.
- **Same-sequence fragment misbinding** — Codex's `part` preservation made a
  v1.0.13 defect visible: `mod(mdg4)`'s first fragment was written at
  another fragment's coordinates, its own locus without a gene row. Reach
  measured first with the rule's own predicate on all six Liftoff inputs: one
  family, every child with exactly one containing fragment.
- The miniprot-candidate guard, the GTF temp-dir leak.
- Two qualification-tool defects: the ladder counted the summary line
  `Errors : 0` as an error (and validated with the environment's installed
  v1.0.13), and dependency evidence crashed on a legacy egg-info install,
  failing 8 tests on the Python 3.10 environment only.

**Admitted, twice.** Cycle 4's CDS fix attached the *unclipped* CDS to one
exon, so containment normalisation later widened that exon across the intron;
its test stopped before serialization. Codex caught it, and the replacement's
tests assert on the written, validated, translated GFF3. And this round's
first cut of the re-enumeration fix (`18351df`) pre-scanned every lifted root,
breaking the dispatcher's lazy, window-bounded root scan; the targeted tests I
ran passed and the full suite failed it on all three Pythons (1,004 roots read
where the window allows 4). `3265d87` tracks yielded ids instead.

## Qualification evidence (`e688ff5`, 2026-09-23)

*The 2026-09-24 requalification of `3ab51cf` is in the section above.*

Every run below imported LiftOn from a detached worktree at the then-frozen
commit and asserted so on load (`build: …/wt_e688ff5_*/lifton/__init__.py e688ff5
dirty=0`); lifts shared cached `-L`/`-M` from the v1.0.12 release validation
so the build is the only difference; long jobs ran in an isolated tmux server.
The shipped code (`lifton/`, `setup.py`, `pyproject.toml`, `MANIFEST.in`,
`lifton.yml`) is unchanged between `e688ff5` and the pushed head — later
commits touch only `notes/` and the two changelogs, neither of which ships.

### Tests and static gates (`e688ff5`)

| gate | result |
|---|---|
| full suite, Python 3.10.21 | **2,485 passed**, 14 skipped, 0 failed (50:37) |
| full suite, Python 3.11.15 | **2,485 passed**, 14 skipped, 0 failed (50:44) |
| full suite, Python 3.12.14 | **2,485 passed**, 14 skipped, 0 failed (50:36) — the `18351df` run's 2,444 (incl. its two lazy-root-scan failures, fixed in `3265d87`) plus 41 new `transl_except` tests; the same 14 skips |
| 24-cell byte-identity matrix + integration (`make test-fast`) | 30 passed, no golden edit |
| fatal flake8 (CI's `E9,F63,F7,F82` over `.`) | 0 |
| Sphinx docs build (pinned Sphinx 9.1.0) | succeeds; 77 warnings, **0 new** vs v1.0.13 (which had 80) |
| `make benchmark-gate` (isolated export of `e688ff5` with its own copy of `work/human_mane`, PYTHONPATH-pinned; the lift's manifest carries the new `transl_except` counters) | **GATE PASS**: 24-cell + integration pytest pass; human_mane protein identity 0.99425 → 0.99542 (the evaluator now reads Sec through for every tool — Liftoff's own score rose too), completeness 0.99764 unchanged, wall 11.7 s vs the baseline's 21.2 s (an old baseline; not a speed claim) |
| CI on the pushed head (`51ebcf7`; later commits are this note's own) | **green**: "Run tests" run 35951581433 — build 3.10 / 3.11 / 3.12 and optional-native 3.10 / 3.11 / 3.12; "Qualify distribution artifacts" run 35951581489 — build, native artifact, compiler-free wheel and sdist on 3.10 / 3.11 / 3.12 / 3.14 |

### The new validator rule, measured before shipping

`cds_exon_containment` fired **0** times on 13 files when first measured:
seven whole-genome LiftOn outputs of the previous build and six RefSeq
reference annotations (rice IRGSP, drosophila, bee, dog, human
`NCBI_RefSeq_no_rRNA`, human RS_2025_08 primary). It fires 0 times on all
seven `e688ff5` whole-genome outputs too (below: 0 errors of any kind).
The references fail other checks (duplicate IDs, organellar trans-splicing —
622 to 26,944 errors each), so the rule is not being starved of odd input;
it simply has no false positives there.

### Whole genome: `3116551` (cycle-4 head) → `e688ff5`

Two layers, each measured on its own: `3116551` → `18351df` is this
round's Step-7 work (the previous qualification), and `18351df` → `e688ff5`
touches only transcripts whose reference declares `transl_except` (next
section: 0 rows outside them on every genome).

| genome | rows only in `3116551` → only in `e688ff5` | genes / transcripts | identity worse / better | status |
|---|---|---|---:|---|
| bee | 502 → 581, all in the 46 declaring transcripts | 12,390 / 28,080, unchanged | 0 / 34 | success |
| rice | 688 → 813: 5 rows `nad5` rebound to its own sequence, `rps12`/`nad1`/`nad5` source `part`; the rest in 109 declaring transcripts | 34,374 / 55,002, unchanged | 0 / 84 | success |
| drosophila | 4,008 → 4,072: 35 rows `mod(mdg4)` (5 fragments at their own spans, 31 transcripts under the fragment that contains them); the rest in 513 declaring transcripts | 16,585 / 33,623, unchanged | 0 / 482 | success |
| human → zebrafish | 325 → 478, all in declaring transcripts and one GPX1 isoform | 15,179 → **15,186** / 69,012 → **69,030** (seven selenoprotein genes recovered) | 0 / 18 | success |
| dog → cat | 1,613 → 1,960: `gene-LOC102156326` restored (+37 transcripts); the rest in 247 declaring transcripts | 32,396 → 32,397 / 90,191 → 90,228 | 1 (LCE6A, −0.007) / 43 | `partial_success`, 396 failures → **success, 0** |

"Transcripts lost/added" in the raw comparison (bee 6/6, rice 17/17,
drosophila 4/4, dog → cat 2/39, human → zebrafish 1/19) are keyed by ID *and*
span: the paired ones are the same IDs whose span moved when a declared stop
was read through; no transcript ID disappeared.

`gff3-validate` (v1.0.14, uncapped) on every `e688ff5` output: **valid, 0
errors** on all five (bee 528,245 rows, rice 665,442, drosophila 395,268,
human → zebrafish 1,730,264, dog → cat 1,815,853), and on the
human → zebrafish `--no-rescue-second-locus` arm (1,714,687 rows). Warning classes are
identical to `18351df`; only the `non_cds_phase` count moves, within the
declaring transcripts: bee 622 → 553, rice 615 → 535, drosophila 1,711 →
1,697, dog → cat 16,072 → 16,047, human → zebrafish 795,401 → 795,465 (the
exons of the rescued selenoproteins). Run status `success`, 0 failures on all
five.

**Threading contract** (`-t 1`, the default, vs `-t 8`), `e688ff5`: rice
identical (187,153,302 B), drosophila identical (162,660,686 B; exercises the
parent overlay's threaded reopen), dog → cat identical (523,529,843 B; where both of this
round's Step-7 fixes act).

### `transl_except`

The gate on each pair of outputs (`v14_q3/transl_except_gate.py`): every
changed row must belong to a transcript whose *reference* declares
`transl_except`; every emitted `transl_except` must lie inside its own CDS,
in frame, with a stop codon in the target genome for Sec/Pyl/Other.

| cell | changed rows outside declaring transcripts | Sec identity before → after | worse | codon check |
|---|---:|---|---:|---|
| MANE chr22 (GRCh38 → CHM13), `3265d87` → `e688ff5` | **0** of 91 | 0.834 → **0.998** (3 of 3 better; SELENOM 0.514 → 1.000) | 0 | 4 of 4 ok |
| CHM13 whole genome (RS_2025_08), `18351df` → `e688ff5` | **0** of 1,671 | **0.666 → 0.998** (53 of 53 better; below 0.9: 36 → 0; lowest 0.980) | 0 | 127 of 127 ok |
| five whole-genome pairs (vs the `18351df` arms) | **0** | see below | 1 (`Other`, −0.007) | 881 of 881 ok |

Whole-genome arms, `18351df` → `e688ff5`, `-t 8`, same cached `-L`/`-M`
(`v14_q3/te/run_all.sh`; identity is over models present in both builds):

| genome | changed rows (transcripts) | outside | identity, declared kinds | codons |
|---|---|---:|---|---|
| bee | 1,083 (46) | 0 | `Other` 0.642 → 0.960 (37 of 40 better); Sec 0.79 → 1.00; Cys/Tyr/Lys 2 of 2 better | 62 / 62 |
| rice | 1,491 (109) | 0 | `Other` 0.660 → 0.987 (101 of 107 better) | 131 / 131 |
| drosophila | 8,010 (513) | 0 | stop-readthrough isoforms (`Other`) 0.852 → 0.965 (482 of 486 better); Sec 0.758 → 0.906 (4 of 4) | 563 / 563 |
| human → zebrafish | 803 (45) | 0 ¹ | Sec 0.543 → 0.650 (12 of 12 better) **+ 18 Sec transcripts that were missing** (mean 0.613) | 45 / 45 |
| dog → cat | 3,340 (247) | 0 | Sec 0.575 → 0.859 (23 of 27 better); `Other` 20 better, **1 worse** ² | 80 / 80 |

No model of any kind got worse except the one in note 2; Met, TERM and the
amino-acid-over-stop kinds are unchanged wherever they did not improve.

¹ Plus 8 rows of GPX1's `NM_001329455.2`, an isoform declaring nothing: before
the fix the rescue placed GPX1 through that isoform at zebrafish `gpx1b`,
because its Sec-declaring transcripts scored below the rescue floor; it now
places GPX1 through `NM_000581.4` at `gpx1a` and the isoform follows its gene.
Both are orthologs of human GPX1.

**The 18 new transcripts are real.** They are seven genes the lift used to
miss entirely — GPX4, DIO1, DIO2, DIO3, SEPHS2, SELENOT, SELENOM; genes
15,179 → 15,186, exactly these seven — plus new isoforms of GPX1 and
SELENOF. They come from the miniprot-only rescue,
whose protein-identity floor a selenoprotein could not clear while identity
stopped at the Sec codon. Checked against zebrafish's own GRCz11 annotation,
every one lies on the corresponding zebrafish gene: GPX1 → `gpx1a`, GPX4 →
`gpx4b`, DIO1 → `dio1`, DIO2 → `dio2`, DIO3 → `dio3b`, SEPHS2 → `sephs3`
(zebrafish has no `sephs2`; `sephs1` is the non-Sec paralog), SELENOF →
`selenof`, SELENOT → `selenot1a`, SELENOM → `selenom`. Existing models improved
the same way (GPX2 0.147 → 0.723, GPX3 0.093 → 0.525).

² LCE6A `XM_022404378.2`, 0.407 → 0.400: RefSeq's dog model reads through a
genomic stop (`aa:Other`). Before, the model was cut there and the ORF rescue
replaced it with an ORF ending at that stop (0.407). Now the lifted model
reads through and scores 0.400, and the ORF rescue replaces a CDS only for a
gain above 1 % (`__find_orfs`, `threshold_orf`, unchanged), so the reference's
two-exon structure is kept, with its `transl_except` at the target codon.

Rice's chloroplast `rpl2` (ACG start, CDS under the gene) lost its
`transl_except` instead of having it rewritten: LiftOn cannot align that model
to a reference protein (`mutation=no_protein`, unchanged since `3116551`), so
the Met cannot be placed — and the value it carried before was a
chloroplast (NC_001320.1) coordinate on a CDS lifted to CP132244.1.

Hermetic fixture (tests/test_transl_except.py) on the previous build
reproduces the report exactly — the SEPHS2-shaped gene's CDS start moved and
scored 0.885, a split-codon selenoprotein scored 0.519, both carrying
reference-coordinate `transl_except` — and passes on `e688ff5`.

### Second-locus rescue (default on)

| check | result |
|---|---|
| human → zebrafish vs zebrafish's own GRCz11 RefSeq | **655 of 689** placements match a distinct protein-coding target gene (≥50 % reciprocal CDS overlap); shifted-locus null over 1,000 replicates: mean 0.57, max 5; **p = 0.001**; 0 rows of the rescue-off output lost. Rerun on the `e688ff5` outputs (sha256 `09b85aca…` on, `aec977e0…` off): unchanged from the `2bae0b9` result, although 18 Sec transcripts moved or appeared — none is a second-locus placement. |
| rice → sorghum vs sorghum's own RefSeq (Codex session) | 51 of 53, null max 2 of 1,000, p = 0.001, 0 rows lost |
| eight-cell safety ladder, `e688ff5` | **8/8 pass** — 0 lost, 0 regressed, 0 overlapping, validity 0 → 0 in every cell; placements identical to the 09-19 promotion run |

### A gain this release already carries, attributed

On the ladder's rice → sorghum subset the rescue-off arm itself moved between
the 09-19 build and today: **1,592 transcripts gained protein identity (mean
+0.119 among those that changed), 15 lost a little (largest −0.064), none
dropped.** Rerunning that arm on the builds either side of cycle 3's P2
(`1963d5b`, miniprot's redundant `stop_codon` no longer ingested as an exon)
reproduces the 09-19 output byte-for-byte before and today's output
byte-for-byte after: the whole effect is that one commit. P2's own A/B
(`notes/overlapping_exons_2026-09.md`) covered five other genomes and did not
include this cross-species pair, so this size of effect was not on record.

### CHM13

Regenerated on `e688ff5` with the same recipe as the staged file (cached
Liftoff, fresh miniprot, `-t 16`), 1 h 09 m, peak RSS 30.2 GiB, status
`success`, 0 failures: `/ccb/salz3/kh.chao/lifton_chm13_v1014_final/`
(`REPORT.md` there). **Staged, not published.**

- The file staged until now (`lifton_chm13_regen2`, cycle 3) is
  byte-identical to the `18351df` regeneration, so the whole difference is
  the `transl_except` fix: 1,671 rows, all in the 99 transcripts that declare
  it; gene and transcript ID sets identical (42,689 genes, 131,823 mRNA).
- Selenoproteins 0.666 → 0.998 (53 of 53, all 25 genes; none below 0.9);
  `Other` 0.934 → 0.999; Ser/Trp 3 of 3 to ≥ 0.99; Met, TERM unchanged;
  nothing worse.
- SEPHS2: Liftoff's full model, identity 1.000, CDS chr16:30,830,620–30,831,966,
  `transl_except=(pos:complement(30831787..30831789),aa:Sec)` — a TGA in
  CHM13. The reporter's file carried GRCh38's `30445548..30445550`.
- 127 values written, 127 placed correctly; 2 not written by design
  (LOC102724117: no stop there in CHM13; MUC19: a miniprot model the codon
  cannot be placed on).
- `gff3-validate`: valid, 0 errors (12,031 warnings; only `non_cds_phase`
  moved, 1,262 → 1,254). Issue #26 overlaps 0 / 0, issue #16 duplicate exon
  IDs 0 (NOC2L 19 of 19 distinct), mitochondrial CDS 13 of 13.

### Packaging (`e688ff5`)

| check | result |
|---|---|
| sdist / wheel | built from `git archive e688ff5`: wheel `lifton-1.0.14-py3-none-any.whl` sha256 `8e518c3c…c2979`, sdist `lifton-1.0.14.tar.gz` sha256 `f700460b…b37e` |
| `twine check --strict` | pass (both) |
| wheel `lifton/` vs the frozen tree | 95 files, **0 differ**; the 12 not shipped are vendored Liftoff's own tests |
| install with `CC=/bin/false` (wheel on 3.10 / 3.11 / 3.12, sdist on 3.11) | all four succeed with no compiler; the only sdist-built dependency is pure-Python `interlap` (3.11/3.12 build its wheel; 3.10's older pip uses the legacy `setup.py install`), and the sdist install also builds `lifton` itself; `mappy` absent as intended; `lifton -V` = v1.0.14 |
| installed console script, fresh chr22 lift (minimap2 2.28-r1209, miniprot 0.13-r248, `-copies`) | 4 of 4 exit 0, `gff3-validate` exit 0; 77,854 rows, 894 genes, 2,801 mRNA; **all four outputs byte-identical** (md5 `4d98f5c8…`). vs the `18351df` smoke lift: 243 rows differ, all in the 9 chr22 transcripts that declare `transl_except` (Sec 0.914 → 0.998, 6 of 6 better; 9 of 9 codons placed) |
| dependencies vs v1.0.13 | unchanged (so the Bioconda bump is version + sdist hash) |

## Known limits — shipped as they are, on purpose

- **Degenerate lifts are emitted as Liftoff made them.** The dog → cat gene
  this release recovers has five mRNAs whose lifted CDS is a single 3-bp
  fragment; they carry `status=no_ref_protein`, the existing label for "no
  protein to compare" (a legacy name — it also covers an empty *lifted*
  protein). Emitting them is faithful to Liftoff; dropping 37 transcripts for
  it was not.
- **`genes_emitted_without_children` counts rebound trans-spliced fragments.**
  A fragment whose transcripts moved to the fragment that contains them is now
  emitted without children, and the counter reads that as a loss: drosophila
  0 → 1 (`mod(mdg4)`), rice 2 → 3 (`nad5`). No transcript was lost in either;
  the counter cannot tell the difference yet.
- **Exons of miniprot-derived models carry a phase.** They are cloned from
  miniprot's CDS rows and keep its phase, score and alignment attributes, so
  `gff3-validate` warns `non_cds_phase` (GFF3 practice is `.` off CDS): on
  human → zebrafish 795,401 of 828,497 exon rows — 762,075 from the
  miniprot-only rescue, 29,635 from Step 8, 3,691 from candidate 3 — and
  16,072 on dog → cat, 622 on bee. A warning, not an error: no translation
  reads an exon's phase, and the count is identical in the cycle-4 output.
  Deferred rather than fixed at the end of qualification because clearing it
  rewrites most exon rows of every cross-species lift; the next version should
  do it as its own change, gated on every column but exon column 8 being
  byte-identical.
- **`transl_except`, where it stops.**
  - miniprot's input keeps `*` at a selenocysteine (`proteins.fa` is
    unchanged on purpose). In every run measured miniprot aligned through it —
    all 30 Sec models on human → zebrafish, 18 of them from the rescue, carry
    a placed `transl_except` — but nothing forces it to; a miniprot model that
    stopped there would be scored correctly and stay short. Writing `U` into
    `proteins.fa` is untested with miniprot (next version).
  - ~~The ORF scan cannot read through a declared codon~~ — fixed in the
    second scan (M3): dog → cat LCE6A is 0.467, above both earlier values.
  - A model LiftOn cannot align to a reference protein (`no_protein`) gets its
    `transl_except` removed, not rewritten (rice chloroplast `rpl2`): a value
    is written only where it can be placed on the model.
  - The benchmark evaluator now reads declared stops through, so its scores
    change for those transcripts; archived results were not re-scored.
- **The CDS-split path is not reached by the corpus.** Its correctness rests
  on tests that serialize, validate and translate the result, on both strands.
- **Reference-keyed recall cannot see the second-locus gain.** The ladder is a
  safety gate; whether placements are real needs the target's own annotation,
  which exists for two cells (human → zebrafish, rice → sorghum), both above.
- **Deferred performance** (measured, next version): gffutils has no batched
  `children()`, so Step 7, Step 8 and the rescue prefetch query row by row
  (~34 % of Step-7 dispatch on rice); `__find_orfs` is 8 % of mammalian
  dispatch (`notes/step7_profile_2026-09-21.md`).
- **77 pre-existing Sphinx warnings** (theme options, heading underlines,
  duplicate section labels in the tutorials). None new; three fewer than
  v1.0.13.
- **Qualification environments.** The Python 3.10 and 3.11 environments
  still live in `/tmp` (created by the Codex session); the 3.12 one lost files
  to tmp cleanup on 2026-09-24 and was rebuilt as
  `/ccb/salz3/kh.chao/lifton_improve/qual_envs/py312` with the 3.11 package
  set. Move the other two there before the next release. Dependency evidence
  for a legacy egg-info install (interlap on 3.10) is recorded as weaker
  (`inventory: egg-info SOURCES.txt`) rather than failing.
- **A gene nested in a neighbour's intron can be suppressed** once that
  neighbour's model is complete. Step 8 and the rescue refuse a miniprot model
  whose span overlaps an emitted gene's span, not its exons — a rule older
  than v1.0.14. v1.0.14 completes more models (`1963d5b`), so it reaches a few
  more nested genes: 12 of the 13 transcripts the 34 subsets give up, 2 of the
  15 on the 8 genomes. Suppressing by exon overlap is a recall feature with its
  own duplicate risk (Iterations 13 and 15); next version.
- **Two CDS rows sharing a base are refused**, counted as
  `cds_spanning_exons`: 1 transcript in the 34 subsets (C. elegans →
  C. briggsae W07G4.3), 0 on the genomes, CHM13 and MANE. Keeping it needs
  the double-counted base resolved, which changes the protein.
- **Six arabidopsis mitochondrial cosRNAs** (15–46 bp `ncRNA` children of
  coding genes) are absent from both v1.0.13's and v1.0.14's arabidopsis →
  rice output, uncounted; v1.0.13 also logged them as failures.
- **Rescue isoform scoring trades time for memory.** `--rescue-max-inflight`
  (default 8,192) restarts the worker pool per batch: human → chicken scores
  its 46,241 isoform jobs 55 s slower than v1.0.13 but peaks 1.1 GiB lower;
  `--rescue-max-inflight 0` restores one batch.

## Release runbook (not executed — each step is an outward action needing sign-off)

State at handoff: `v1014-integration` = the frozen code `3ab51cf` plus
commits touching only `notes/` and the two changelogs (neither ships in the
wheel; the changelogs are checked by `test_packaging_metadata`),
pushed, CI pending (runs on push).
`main` = `devel` = `b2fe59f` (v1.0.13); the branch is 77 commits
ahead and 0 behind, so both merges are fast-forwards.

1. **Date the release.** Replace the provisional `2026-09-22` in all three
   places `tests/test_packaging_metadata.py` requires to agree:
   `CITATION.cff` (`date-released`), `CHANGELOG.md` (`## [1.0.14] - DATE`),
   `docs/source/content/changelog.rst` (`v1.0.14 (DATE)`). Run
   `pytest tests/test_packaging_metadata.py`, commit
   `chore(release): date v1.0.14 as DATE` (the v1.0.13 precedent is `b2fe59f`),
   push, wait for CI.
2. **Fast-forward the branches.**
   ```
   git checkout devel && git merge --ff-only v1014-integration && git push origin devel
   git checkout main  && git merge --ff-only devel            && git push origin main
   ```
   Pushing `main` rebuilds the docs site (`.github/workflows/docs.yml` →
   khchao.com/LiftOn). Check it rendered v1.0.14.
3. **Optional dry run.** Run `publish.yml` by hand (`workflow_dispatch`): it
   builds and qualifies the exact sdist + wheel (`packaging.yml`) and publishes
   to **TestPyPI** only. Install from TestPyPI in a clean venv and lift chr22.
4. **Tag.** `git tag -a v1.0.14 -m "LiftOn v1.0.14" && git push origin v1.0.14`
5. **GitHub Release.** `gh release create v1.0.14 --title "LiftOn v1.0.14"
   --notes-file notes/release_notes_v1.0.14.md`. Publishing it triggers
   `.github/workflows/publish.yml` (`release: published`), which builds and
   qualifies the sdist + wheel via `packaging.yml` and uploads them to PyPI
   (OIDC trusted publishing).
6. **Verify all four surfaces** (the v1.0.10 lesson): tag, Release, PyPI
   (`pip download lifton==1.0.14 --no-deps --no-binary :all: --no-cache-dir`),
   and `main` ancestry (`git merge-base --is-ancestor v1.0.14 origin/main`).
   Then a clean-venv `pip install --no-cache-dir lifton==1.0.14` and one real
   lift (the chr22 example) from the installed package.
7. **Bioconda PR #66594** (`Kuanhao-Chao/bioconda-recipes`, branch
   `add-lifton-1.0.9`, currently "Add lifton 1.0.13", open, awaiting review).
   Runtime dependencies are unchanged since 1.0.13, so the recipe change is
   `version: 1.0.14` + the `sha256` of the **PyPI** sdist from step 6 (not
   the qualification build's hash — CI rebuilds from the tag). Retitle the PR
   "Add lifton 1.0.14". It is still unmerged, so bumping it avoids a
   1.0.13 → 1.0.14 autobump that would ship 1.0.13's defects first.
   Reviewer nit (2026-09-20): the Sep 18 comment's second `#78` links to
   bioconda-recipes #78; edit it to
   `https://github.com/Kuanhao-Chao/LiftOn/issues/78`.
8. **CHM13 annotation.** The v1.0.14 regeneration is staged at
   `/ccb/salz3/kh.chao/lifton_chm13_v1014_final/` (`REPORT.md` there; see the
   CHM13 section). `../lifton_chm13_v1014/` and `../lifton_chm13_regen2/` are
   superseded. Publishing
   = copy to `/ccb/salz7-data/ftp.ccb/pub/data/LiftOn/` as
   `JHU_LiftOn_v1.0.14_chm13v2.0.gff3` (+ a statistics sheet like the
   v1.0.12 one), replace `human_refseq/lifton.gff3` (currently the v1.0.12
   file), and point the README/docs link at it (a docs commit; can ride
   step 1). The derived tracks under `lifton_chm13_2026/` (`lifton.bb`,
   `mutations/`, `visualization/`) are older still.
9. **Issue #16** ("Provided CHM13 file has incorrect exons", open). Its reply
   (2026-07-30) says the tool bug is fixed and the posted file is not yet
   replaced. After step 8, reply that the posted file has been replaced, with
   the NOC2L check (19 exons, 19 distinct IDs; 0 duplicate exon IDs in the
   file) and the #26 check (0 overlapping exons or CDS). Use `gh api repos/Kuanhao-Chao/LiftOn/issues/16/comments -f body=...`
   (`gh issue comment` is broken by the Projects-classic deprecation).
10. **Reply to the `transl_except` reporter** (email, not a GitHub issue):
    `notes/reply_simon_transl_except.md`. Fill its availability line once the
    release exists (steps 5–6), then send it yourself. It is about their own
    MANE → CHM13 lift, so it does not depend on step 8.

## Appendix A — Release readiness — v1.0.14 (cycle 3), as committed in `4d532b2`

Written 2026-09-21 on branch `v1014-integration`. Cycles 1 and 2 landed 23
commits; this cycle adds three. Two of them are correctness fixes for defects
that shipped in every release; the third is a measured speed win on the
windowed aligner.

Nothing here has been pushed, tagged, released or published.

### What changed

| | what | class |
|---|---|---|
| P1 | `lifton_add_trans_exon_cds` asks for **level-1** exons, so `--threads 1` matches `--threads N` | output-corrective |
| P2 | no transcript emits overlapping exons: miniprot's redundant `stop_codon` is no longer ingested as a second exon, and `update_cds_list` reconciles a rebuilt exon against the one it ran into | output-corrective |
| P3 | `windowed_align._unique_anchors` indexes the reference only over the query's k-mers | byte-neutral |

Details and the reasoning behind each are in
`notes/threading_exon_divergence_2026-09.md`,
`notes/overlapping_exons_2026-09.md` and
`notes/windowed_anchor_construction_2026-09.md`.

### What was wrong the first time

Two claims in this cycle's own working notes had to be withdrawn after
measurement, and both are worth keeping visible.

**The P1 blast radius was counted in the wrong database.** The audit located
the divergence correctly and then counted how often it could fire by scanning
the *reference* annotation (1,915 human loci, 46 on chr22). The query reads
Liftoff's *output*, which contains none of that shape, because Liftoff does not
lift the nested miRNA. The real trigger is a different one entirely — RefSeq's
organellar convention — and it appears on rice (17) and arabidopsis (7), not on
human at all. A whole-genome human → CHM13 A/B was run on the strength of the
wrong count and correctly showed no change.

**The P2 gate was specified too strictly.** "Protein identity unchanged on
every already-valid transcript" failed on the first run — two rice transcripts
moved. Both had moved *up*, because repairing the terminal CDS let the miniprot
candidate be scored on its real sequence. The invariant that matters is that no
already-valid transcript gets worse, and the gate was corrected to that rather
than the result being explained away.

### Verification

Every A/B arm ran in a detached tmux session from a build pinned to an explicit
worktree, and asserted on load which `lifton/__init__.py` it had imported.
Paired arms shared one cached `-L`/`-M` so the build is the only difference.
Where an arm was started before a pinned worktree existed, it was re-run from
the pinned build and the two outputs compared byte-for-byte (identical).

#### P1 — does `--threads 1` equal `--threads N`?

| | before: `-t 1` vs `-t 8` | after: `-t 1` vs `-t 8` |
|---|---|---|
| rice | **DIFFER** — 187,083,978 vs 187,068,200 bytes | identical |
| human → CHM13 | identical (1,237,342,491 both) | identical |

On rice the pre-fix **serial** arm is the only one of the four that differs:
`before -t 8`, `after -t 1` and `after -t 8` are byte-identical to each other.
Seventeen transcripts lose 35 exon rows and 35 CDS rows; all 17 carried
duplicate-coordinate exons before and none do after. `-t 1` is the default, so
the default path was the wrong arm.

Human → CHM13 is inert, which the Liftoff-output scan predicted (0 features
carrying the shape) and three measurements confirm: pre-fix `-t 1` equals
pre-fix `-t 8`; the pinned P1 build at `-t 8` is byte-identical to the pinned
pre-P1 build at `-t 8`; and the P1 build at `-t 1` closes the pair.

#### P2 — overlapping exons, five whole genomes

Both arms pinned, one shared cached `-L`/`-M` per pair.

| | human → zebrafish | drosophila | CHM13 | rice | bee |
|---|---:|---:|---:|---:|---:|
| overlapping-exon transcripts | 96 → **0** | 47 → **0** | 27 → **0** | 13 → **0** | 13 → **0** |
| overlapping-CDS transcripts | 95 → **0** | 16 → **0** | 13 → **0** | 4 → **0** | 8 → **0** |
| genes / transcripts | unchanged | unchanged | unchanged | unchanged | unchanged |
| already-valid transcripts worse | 0 | 0 | 0 | 0 | 0 |
| already-valid transcripts better | 4 | 28 | 4 | 2 | 3 |
| pairs crossing a strand/seqid | 0 | 0 | 0 | 0 | 0 |
| "spans 2 exons" warnings | — | 2,632 → 0 | 8,546 → 0 | 4,468 → 0 | 3,614 → 0 |

(rice's 13 is what remains after P1 removed the 17 it was responsible for.)

#### P3 — is the aligner change output-safe?

Whole-genome dog → cat, both arms pinned to frozen worktrees differing only by
this change: **523,466,820 bytes, byte-identical.**

#### Suite

2,399 passed, 2 skipped, 0 failed (2,358 at the start of the cycle). 24-cell
matrix green with no golden edit. Fatal flake8 clean.

#### P5 — is the drop ledger visible on real data?

The counter reaches `run_manifest.json` on every real run, with all seven
classes recorded including the new `hierarchy_depth_exceeded`, and the
end-of-run summary correctly stays silent when nothing was dropped.

It has still not been seen firing outside a test, because on these corpora
nothing is dropped. The 550 `Skipping … was not found` lines the plan expected
from rice were the `-copies` resolution bug, fixed in v1.0.12. That is the
right answer for these inputs, not a gap in the instrument — but it does mean
the classes are exercised only by unit tests.

### Known, not fixed

* **`nad5` in rice** — `rna-OrsajM_p05` is written on `CP132246.1` while its
  gene `gene-OrsajM_p05` stays on `CP132245.1`. A gene and its transcript on
  different sequences is invalid GFF3. It is a trans-spliced mitochondrial
  model (`exception=trans-splicing`), present identically before and after this
  cycle's changes. Surfaced by the new validator check; not caused by it and
  not fixed by it.
* **Two disjoint coding blocks under one overlapping exon pair** — an exon
  holds one CDS, so `reconcile_overlapping_exons` refuses rather than inventing
  coding sequence. No such pair occurs in the five genomes measured.

### A process failure worth recording

Mid-cycle I rewrote `p2_ab/arm.sh` in place to add a fifth genome. Two
`human → zebrafish` arms were still running, and bash reads a script
incrementally from an open file descriptor: `open(path, 'w')` truncates and
rewrites the **same inode**, so every byte offset after the insertion shifted
and both shells resumed at a misaligned position, re-executing the tail of the
script and starting a second lift on top of a finished one.

The artifacts survived, and were verified rather than assumed:

* `out.gff3` was byte-unchanged from the copy taken the moment the problem was
  noticed, so the second run never reached the publish step;
* there is no `*.partial.gff3` and `run_manifest.json` records
  `status: success` for both arms — `OutputTransaction` publishes only on
  success;
* the `before` arm is **byte-identical to `c1_ab/on`**, an independent,
  earlier, complete run of the same configuration.

That third check also independently confirms the P1 analysis: `c1_ab/on`
predates P1 and `before` includes it, and they are the same file, which is what
the Liftoff-output scan predicted for a genome carrying none of the affected
shape.

The rule this cycle already had — *pin every A/B arm to an explicit build* —
did not cover the driver script itself. A running script is as much live state
as a running tree. The arm scripts are now read-only, and a new variant goes in
a new file.

**A second one, same family.** The human → CHM13 A/B began as a four-arm serial
driver. To parallelise it I wrote empty placeholder `out.gff3` files so the
driver's `[ -f out.gff3 ]` guard would skip the arms I was moving, and gave the
parallel arms an `rm -f out.gff3` so they would not inherit a placeholder. The
two interact: the parallel arm deleted the placeholder, the driver reached that
arm before the parallel one had published, saw no file, and ran it again — from
the **live** tree, which by then carried P2 and P3. Its `after_t1` output
therefore had 0 overlapping-exon transcripts where its `after_t8` sibling had
27, and the difference read at first like a second threading divergence.

Diagnosis came from the driver's own log (`[after_t1] build: .../src/...`,
where the pinned arms log a worktree path) and the output timestamp, two hours
after the parallel arm had finished.

**The claim was then rebuilt on arms that are actually pinned, and it holds.**
Commit `3a92400`'s human → CHM13 row rests on these three measurements, none of
which involve the discarded arm:

| | |
|---|---|
| pre-fix `-t 1` vs pre-fix `-t 8` | identical, 1,237,342,491 bytes (both pinned to the pre-P1 worktree) |
| P1 build `-t 8` vs pre-P1 build `-t 8` | byte-identical — the fix changes nothing here |
| P1 build `-t 1` vs P1 build `-t 8` | byte-identical (`p1_chm13_clean/after_t1` vs `p2_ab/chm13/before`) |

All four arms are the same bytes, which is what the Liftoff-output scan
predicted for a genome carrying none of the affected shape. Anyone reproducing
this from the artifacts on disk should use `p1_chm13_clean/after_t1`, not
`p1_ab/after_t1` — the latter is the discarded arm and is left in place only so
this note can point at it.

The general shape, for the third time this cycle: a guard is only a guard if
nothing else is allowed to change what it tests.

### The regenerated CHM13 annotation

`/ccb/salz3/kh.chao/lifton_chm13_regen2/` — 1,237,332,649 bytes, 1 h 06 m,
peak RSS 30.3 GiB. Same recipe as the cycle-2 regeneration so the two compare
directly. **Staged, not published.**

| | staged (cycle 2) | regenerated (cycle 3) |
|---|---:|---:|
| genes / transcripts | 42,689 / 184,596 | unchanged |
| transcripts with overlapping exons | 27 | **0** |
| transcripts with overlapping CDS | 13 | **0** |
| `"reference model is malformed"` warnings | 8,537 | **0** |
| `gff3-validate` | `False`, 40 errors | **`True`, 0 errors** |

The 40 errors in the previous file were *all* of the kind this cycle fixed (27
`exon_overlap`, 13 `cds_overlap`), so the regenerated annotation is completely
clean. The cycle-1 genetic-code fix still holds: all 13 mitochondrial CDS match
reference length exactly.

This is what lets issue #26 be answered as fixed on both halves rather than
half-fixed — with the caveat, stated in every draft reply, that the **posted**
annotation has not been replaced.
