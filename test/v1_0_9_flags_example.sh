#!/usr/bin/env bash
#
# LiftOn v1.0.9 — user-facing flags example
# ------------------------------------------
# This script demonstrates the LiftOn command-line flags that are most relevant
# in v1.0.9. It is for reference only — DO NOT run it as-is: replace the
# placeholder input paths (ref.gff3 / ref.fa / target.fa) with your own files.
#
# Inputs (same for every invocation below):
#   ref.gff3   : reference annotation (GFF3)   [-g / --reference-annotation]
#   ref.fa     : reference genome   (FASTA)
#   target.fa  : target genome      (FASTA)
# Output:
#   out.gff3   : the lifted annotation for the target genome [-o / --output]
#
# Usage convention:  lifton -g <ref.gff3> <ref.fa> <target.fa> -o <out.gff3>

set -euo pipefail

# ---------------------------------------------------------------------------
# 1) Basic invocation (defaults)
#    Lifts ALL auto-detected gene-like parent types (genes, pseudogenes,
#    ncRNA_genes, ...) AND runs the default-on miniprot-only rescue pass that
#    recovers reference coding genes the DNA lift missed entirely.
# ---------------------------------------------------------------------------
lifton -g ref.gff3 ref.fa target.fa -o out.gff3

# ---------------------------------------------------------------------------
# 2) --gene-only
#    Restore the pre-Iteration-12 GENE-ONLY lift: process only the `gene`
#    hierarchy and drop every other gene-like top-level type (pseudogenes,
#    ncRNA_genes, structured mobile elements). Use this for manuscript
#    reproduction or when you want strictly the `gene` partition.
#    (Ignored if you pass -f/--features explicitly — your file always wins.)
# ---------------------------------------------------------------------------
lifton -g ref.gff3 ref.fa target.fa -o out.gene_only.gff3 \
    --gene-only

# ---------------------------------------------------------------------------
# 3) --no-miniprot-rescue
#    Opt OUT of the default-on miniprot-only rescue pass. With this flag LiftOn
#    will NOT emit a miniprot-only model for a reference coding gene the DNA lift
#    missed entirely, reproducing the pre-v1.0.9 (pre-Iteration-23) behaviour.
#    Use it to keep output comparable to v1.0.8 / earlier baselines.
#    (The reverse, --miniprot-rescue, is a no-op alias since rescue is now the
#    default.)
# ---------------------------------------------------------------------------
lifton -g ref.gff3 ref.fa target.fa -o out.no_rescue.gff3 \
    --no-miniprot-rescue

# ---------------------------------------------------------------------------
# 4) --validate-output  (optionally with --validate-verbose)
#    After writing out.gff3, re-validate the file for GFF3 format correctness
#    and feature-hierarchy/CDS-phase/containment sanity, printing a structured
#    report to stderr. Add --validate-verbose to also print warnings (not just
#    errors). This affects the exit code, not the output bytes.
# ---------------------------------------------------------------------------
lifton -g ref.gff3 ref.fa target.fa -o out.validated.gff3 \
    --validate-output --validate-verbose

# ---------------------------------------------------------------------------
# 5) --threads 8 --locus-pipeline
#    Speed up the per-locus Step 7 work by fanning it out across 8 threads via a
#    ThreadPoolExecutor (the --locus-pipeline scheduler). Output is emitted in
#    submission order, so --threads N is BYTE-IDENTICAL to --threads 1 — this
#    changes scheduling/wall-clock only, not the algorithm or the output bytes.
#    Both flags are required together. Works on the default backend (no --native
#    needed). miniprot's own -t also scales with --threads.
# ---------------------------------------------------------------------------
lifton -g ref.gff3 ref.fa target.fa -o out.parallel.gff3 \
    --threads 8 --locus-pipeline

# ---------------------------------------------------------------------------
# (Optional) Combine several flags in one run, e.g. gene-only + parallel +
# output validation:
# ---------------------------------------------------------------------------
# lifton -g ref.gff3 ref.fa target.fa -o out.combined.gff3 \
#     --gene-only --threads 8 --locus-pipeline --validate-output
