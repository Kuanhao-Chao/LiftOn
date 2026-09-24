# Draft reply — transl_except / selenocysteine (not sent)

*For Kuan-Hao to review and send. Numbers from the v1.0.14 qualification; see
`notes/release_readiness_v1.0.14.md`.*

---

Hi Simon,

Thank you — you were right, and it was worse than one gene.

LiftOn ignored `transl_except`. Both the reference protein and the lifted
protein carry a stop (`*`) at a selenocysteine's UGA, and LiftOn's protein
identity counted matches only up to the first stop in the lifted protein. So
an identical selenoprotein scored as if it ended at its selenocysteine (for
SEPHS2, 60 of 449 residues), the UGA was also called a gained stop codon,
and the correct lifted model lost to a shorter one — for SEPHS2, a miniprot
model missing the first 118 residues. On a GRCh38 → CHM13 RefSeq lift, all
25 human selenoprotein genes were affected: 53 transcripts at a mean protein
identity of 0.666. A second problem: the `transl_except` written to the
output kept the GRCh38 coordinates.

The fix will be in v1.0.14:

- A codon the annotation declares to read through a stop (selenocysteine,
  pyrrolysine, stop readthrough, or an amino acid declared over a stop) is
  read through when LiftOn scores models, calls variants and chains
  Liftoff/miniprot. SEPHS2 now keeps its full lifted CDS
  (chr16:30,830,620–30,831,966 on CHM13, with both UTRs) at protein identity
  1.000, and the 53 selenoprotein transcripts go from a mean of 0.666 to
  0.998, none below 0.98.
- `transl_except` is written in the lifted model's own coordinates — for
  SEPHS2, `transl_except=(pos:complement(30831787..30831789),aa:Sec)`, which
  is the TGA in CHM13 — and left off where the target codon no longer needs
  it (for example a selenocysteine replaced by cysteine).

⟨FOR KUAN-HAO — one line once released: "It's in v1.0.14, released DATE
(`pip install lifton==1.0.14`)." Until then: "It will be in v1.0.14."⟩

If you can share your MANE v1.5 → CHM13 run once it's out, I'd be glad to
confirm SEPHS2 and the other selenoproteins on your exact setup.

Thanks again for the careful report and the screenshot.

Best,
Kuan-Hao
