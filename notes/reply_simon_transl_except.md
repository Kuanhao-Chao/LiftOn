# Draft reply — transl_except / selenocysteine (not sent)

*For Kuan-Hao to review and send. Numbers from the v1.0.14 qualification and
the re-check on the reporter's own command; see
`notes/release_readiness_v1.0.14.md`.*

---

Hi Simon,

Thank you — you were right, and it was worse than one gene.

LiftOn ignored `transl_except`. Both the reference protein and the lifted
protein carry a stop (`*`) at a selenocysteine's UGA, and LiftOn's protein
identity counted matches only up to the first stop in the lifted protein. So
an identical selenoprotein scored as if it ended at its selenocysteine (for
SEPHS2, 60 of 449 residues), the UGA was also called a gained stop codon, and
the correct lifted model lost to a shorter one. That is exactly your
screenshot: SEPHS2 came out as miniprot's model with its first 118 residues
turned into UTR. And the `transl_except` written to the output kept the GRCh38
coordinates.

The fix will be in v1.0.14. I re-ran your exact command
(`lifton -g MANEv1.5.gff -chroms chrom_mapping.txt -copies -sc 0.9
chm13v2.0.fa hg38.p14.fa`) with it:

- SEPHS2 is now the full lifted model — CDS chr16:30,830,620–30,831,966 with
  both UTRs (30,829,870–30,832,113), protein identity 1.000 — and its
  `transl_except` is in CHM13 coordinates,
  `(pos:complement(30831787..30831789),aa:Sec)`, which is the TGA there.
- All 25 selenoprotein transcripts in MANE v1.5: mean protein identity 0.654
  before, 0.999 after (18 were below 0.9; the lowest is now GPX1 at 0.990,
  whose alanine repeat really is shorter in CHM13). The 9 transcripts
  with a non-AUG start keep their scores and now carry their `transl_except`
  at CHM13 coordinates as well.
- As an independent check, I translated each lifted CDS straight from the
  CHM13 sequence, applied the `transl_except` LiftOn wrote, and compared with
  NCBI's own protein records (where selenocysteine is `U`): 27 of 33 are
  identical. The other six differ by one or two residues because CHM13 and
  GRCh38 genuinely differ there (GPX1's alanine repeat is one of them), or,
  in one case, because NCBI's protein itself differs from the GRCh38 genome.
  With v1.0.13 none of the 33 matched.

Two notes:

- `lifton_output/liftoff/liftoff.gff3` is Liftoff's own intermediate output;
  it still carries the GRCh38 coordinates, as in the line you quoted. The
  LiftOn output file is the one with the rewritten values.
- The same fix covers GENCODE/Ensembl annotations, which mark selenocysteine
  as separate `Selenocysteine` (GTF) or `stop_codon_redefined_as_selenocysteine`
  (GFF3) rows rather than `transl_except`. On GENCODE v49 → CHM13 the 71
  selenoprotein transcripts go from a mean identity of 0.66–0.68 to 0.995–0.997.
  (MANE's Ensembl-format GFF marks no selenocysteine at all, so use the RefSeq
  one, as you did.)

⟨FOR KUAN-HAO — one line once released: "It's in v1.0.14, released DATE
(`pip install lifton==1.0.14`)." Until then: "It will be in v1.0.14."⟩

Thanks again for the careful report and the screenshot.

Best,
Kuan-Hao
