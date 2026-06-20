"""Subset builder: reduce each benchmark to ONE chromosome pair.

Produces, under ``work/<id>/subset/``:
  ref.chrom.fa / .fai      subset reference genome (seqid harmonized to GFF)
  tgt.chrom.fa / .fai      subset target genome (chosen syntenic chromosome[s])
  ref.chrom.gff3           subset reference annotation (single ref chromosome)
  ref.proteins.subset.faa  reference proteins restricted to the subset
  subset.manifest.json     chosen chroms, header maps, counts
  synteny/ref_chrom_to_tgt.paf, target_chroms.txt   (cross-species / AUTO)

All steps are idempotent (guarded by work/<id>/.done/subset.done).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# lifton package (installed in lifton_devel) for protein translation parity.
from Bio.Seq import Seq  # noqa: E402


# ---------------------------------------------------------------------------
# small GFF streaming helpers (no gffutils DB build needed for subsetting)
# ---------------------------------------------------------------------------

def _attrs(col9: str) -> dict:
    """Parse a GFF3 column-9 attribute string into a dict (first value wins)."""
    d = {}
    for kv in col9.rstrip(";").split(";"):
        kv = kv.strip()
        if not kv or "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        d.setdefault(k, v)
    return d


def gff_mrna_counts_per_seqid(gff_path: str) -> Counter:
    """Count mRNA features per seqid by streaming the GFF (cheap)."""
    counts: Counter = Counter()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            if parts[2] == "mRNA":
                counts[parts[0]] += 1
    return counts


def choose_ref_chrom(bench: dict) -> str:
    """Pick the reference chromosome to subset to."""
    rc = bench["ref_chrom"]
    if rc != "AUTO_LARGEST_CODING":
        return rc
    counts = gff_mrna_counts_per_seqid(bench["ref_gff"])
    if not counts:
        raise RuntimeError(f"{bench['id']}: no mRNA features found in {bench['ref_gff']}")
    # Most mRNA-dense seqid. Prefer real chromosomes: drop obvious unplaced
    # scaffolds (heuristic: name containing 'scaffold'/'unplaced'/'random'/'_alt').
    def is_scaffold(s: str) -> bool:
        s2 = s.lower()
        return any(t in s2 for t in ("scaffold", "unplaced", "random", "_alt", "chrun"))
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    for seqid, _n in ranked:
        if not is_scaffold(seqid):
            return seqid
    return ranked[0][0]


def subset_gff(gff_path: str, keep_seqid: str, out_path: Path) -> dict:
    """Write only the records on ``keep_seqid``; preserve gff-version and the
    matching sequence-region directive. Return feature-type counts."""
    counts: Counter = Counter()
    protein_accs: set[str] = set()
    acc_to_mrna: dict[str, str] = {}
    cds_parent_protein: list[tuple[str, str]] = []  # (parent_mrna, protein_acc)
    with open(gff_path) as fh, open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    sr = line.split()
                    if len(sr) > 1 and sr[1] == keep_seqid:
                        out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[0] != keep_seqid:
                continue
            out.write(line if line.endswith("\n") else line + "\n")
            ftype = parts[2]
            counts[ftype] += 1
            if ftype == "CDS":
                a = _attrs(parts[8])
                pacc = a.get("protein_id") or a.get("Name")
                parent = a.get("Parent", "")
                if pacc:
                    protein_accs.add(pacc)
                    if parent:
                        cds_parent_protein.append((parent, pacc))
    for parent, pacc in cds_parent_protein:
        acc_to_mrna.setdefault(pacc, parent)
    return {
        "feature_counts": dict(counts),
        "protein_accs": protein_accs,
        "protein_acc_to_mrna": acc_to_mrna,
    }


# ---------------------------------------------------------------------------
# FASTA helpers (samtools faidx)
# ---------------------------------------------------------------------------

def _faidx_seqids(fa: str, samtools: str) -> list[str]:
    fai = fa + ".fai"
    if not os.path.exists(fai):
        subprocess.run([samtools, "faidx", fa], check=True)
    out = []
    with open(fai) as fh:
        for line in fh:
            out.append(line.split("\t")[0])
    return out


def samtools_extract(genome: str, seqids: list[str], out_fa: Path, samtools: str,
                     rename_to: str | None = None) -> None:
    """Extract ``seqids`` from ``genome`` into ``out_fa`` and index it.
    If ``rename_to`` is given (single seqid), the output header is rewritten."""
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fa, "w") as out:
        proc = subprocess.run([samtools, "faidx", genome, *seqids],
                              stdout=subprocess.PIPE, check=True, text=True)
        text = proc.stdout
        if rename_to is not None and len(seqids) == 1:
            lines = text.splitlines(keepends=True)
            if lines and lines[0].startswith(">"):
                lines[0] = f">{rename_to}\n"
            text = "".join(lines)
        out.write(text)
    subprocess.run([samtools, "faidx", str(out_fa)], check=True)


# ---------------------------------------------------------------------------
# synteny: pick target chromosome(s) for the chosen reference chromosome
# ---------------------------------------------------------------------------

def choose_target_chroms(ref_chrom_fa: Path, tgt_genome: str, out_dir: Path,
                         minimap2: str, threads: int, log=print) -> list[str]:
    """minimap2-align the subset reference chromosome to the full target genome,
    return the dominant syntenic target seqid(s) (>= max(10kb, 10% of matches))."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paf = out_dir / "ref_chrom_to_tgt.paf"

    def _run_minimap(preset: str) -> dict:
        with open(paf, "w") as fh:
            subprocess.run([minimap2, "-x", preset, "-t", str(threads),
                            "--secondary=no", tgt_genome, str(ref_chrom_fa)],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        m: dict[str, int] = defaultdict(int)
        with open(paf) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) < 11:
                    continue
                try:
                    resid = int(p[9])  # residue matches
                except ValueError:
                    continue
                m[p[5]] += resid
        return m

    # asm20 targets ~5-20% divergence; distant pairs (e.g. D.mel -> D.erecta) may
    # yield zero hits, so fall back to the more permissive asm10 before giving up.
    matches = _run_minimap("asm20")
    used_preset = "asm20"
    if not matches:
        log("  [synteny] asm20 found no alignments; retrying with asm10")
        matches = _run_minimap("asm10")
        used_preset = "asm10"
    if not matches:
        raise RuntimeError(f"minimap2 produced no alignments of {ref_chrom_fa} to "
                           f"{tgt_genome} (tried asm20, asm10)")
    total = sum(matches.values())
    thresh = max(10000, int(0.10 * total))
    chosen = sorted([s for s, m in matches.items() if m >= thresh],
                    key=lambda s: -matches[s])
    if not chosen:  # fall back to the single best
        chosen = [max(matches.items(), key=lambda kv: kv[1])[0]]
    (out_dir / "target_chroms.txt").write_text("\n".join(chosen) + "\n")
    log(f"  [synteny] {ref_chrom_fa.name} -> target chrom(s): {chosen} "
        f"({100*matches[chosen[0]]/total:.0f}% top, {used_preset})")
    return chosen


# ---------------------------------------------------------------------------
# reference protein FASTA for the subset
# ---------------------------------------------------------------------------

def _read_fasta(path: str):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def build_subset_proteins(bench: dict, sub: dict, ref_gff_subset: Path,
                          ref_fa_subset: Path, out_faa: Path, log=print) -> int:
    """Write the reference protein FASTA restricted to the subset chromosome.

    - miniprot_target_space == 'transcript' (MANE): translate CDS from the subset
      reference genome, header = mRNA id (rna-NM_...). miniprot Target then maps 1:1.
    - 'protein': filter the provided ``ref_proteins`` FASTA to the protein
      accessions present in the subset GFF (NP_/XP_).
    """
    space = bench.get("miniprot_target_space", "protein")
    n = 0
    if space == "transcript" or not bench.get("ref_proteins"):
        # Translate CDS per mRNA from the subset reference genome (gffutils).
        import gffutils
        import pyfaidx
        from lifton.extract_sequence import get_protein_sequence
        db = _build_db(ref_gff_subset)
        fa = pyfaidx.Fasta(str(ref_fa_subset))
        seen = set()
        with open(out_faa, "w") as out:
            # RefSeq/NCBI annotate transcripts as 'mRNA'; Ensembl/gffread use
            # 'transcript'. Iterate both featuretypes so the transcript-space
            # protein extraction is annotation-source agnostic (a feature
            # carries exactly one type, so RefSeq output is unchanged — it has
            # no 'transcript' rows; the seen-set guards the rare GFF that mixes
            # both). Previously hardcoding 'mRNA' produced 0 proteins on the
            # gffread-Ensembl GFF3 -> empty miniprot -> the pair aborted.
            for ftype in ("mRNA", "transcript"):
                for mrna in db.features_of_type(ftype):
                    if mrna.id in seen:
                        continue
                    cds = list(db.children(mrna, featuretype="CDS", order_by="start"))
                    if not cds:
                        continue
                    prot = get_protein_sequence(mrna, fa, cds)
                    if not prot:
                        continue
                    prot = prot.rstrip("*")
                    if not prot:
                        continue
                    out.write(f">{mrna.id}\n{prot}\n")
                    seen.add(mrna.id)
                    n += 1
        log(f"  [proteins] translated {n} MANE/CDS proteins (transcript space)")
    else:
        allowed = sub["protein_accs"]
        with open(out_faa, "w") as out:
            for name, seq in _read_fasta(bench["ref_proteins"]):
                if name in allowed:
                    out.write(f">{name}\n{seq}\n")
                    n += 1
        log(f"  [proteins] kept {n}/{len(allowed)} subset proteins (protein space)")
    return n


# ---------------------------------------------------------------------------
# gffutils DB helper (inference disabled, force rebuild)
# ---------------------------------------------------------------------------

def _build_db(gff_path: Path):
    import gffutils
    dbfn = str(gff_path) + "_db"
    return gffutils.create_db(
        str(gff_path), dbfn=dbfn, force=True, keep_order=True,
        merge_strategy="create_unique", sort_attribute_values=False,
        disable_infer_genes=True, disable_infer_transcripts=True,
    )


# ---------------------------------------------------------------------------
# top-level
# ---------------------------------------------------------------------------

def build_subset(bench: dict, work_dir: Path, tools: dict, threads: int = 8,
                 force: bool = False, log=print) -> dict:
    sub_dir = work_dir / "subset"
    done = work_dir / ".done" / "subset.done"
    manifest_path = sub_dir / "subset.manifest.json"
    if done.exists() and not force and manifest_path.exists():
        log(f"  [subset] cached ({bench['id']})")
        return json.loads(manifest_path.read_text())

    sub_dir.mkdir(parents=True, exist_ok=True)
    samtools = tools["samtools_bin"]
    minimap2 = tools["minimap2_bin"]

    # --- validate inputs exist, are readable and non-empty (catches stale
    #     0-byte files, mis-set paths, and unreadable target genomes) ---
    required = [("ref_genome", bench["ref_genome"]), ("ref_gff", bench["ref_gff"]),
                ("tgt_genome", bench["tgt_genome"])]
    if bench.get("ref_proteins"):
        required.append(("ref_proteins", bench["ref_proteins"]))
    for label, path in required:
        if not os.path.exists(path):
            raise RuntimeError(f"{bench['id']}: {label} missing: {path}")
        if not os.access(path, os.R_OK):
            raise RuntimeError(f"{bench['id']}: {label} not readable: {path}")
        if os.path.getsize(path) == 0:
            raise RuntimeError(f"{bench['id']}: {label} is empty (0 bytes): {path}")

    ref_chrom = choose_ref_chrom(bench)
    log(f"  [subset] ref chromosome = {ref_chrom}")

    ref_fa = sub_dir / "ref.chrom.fa"
    tgt_fa = sub_dir / "tgt.chrom.fa"
    ref_gff = sub_dir / "ref.chrom.gff3"
    ref_faa = sub_dir / "ref.proteins.subset.faa"

    # --- subset reference GFF (single chromosome) ---
    gff_info = subset_gff(bench["ref_gff"], ref_chrom, ref_gff)
    gff_seqids = {ref_chrom}

    # --- subset reference genome FASTA, harmonize header to the GFF seqid ---
    if bench.get("ready_subset"):
        # Ready chr22 FASTAs already use the GFF seqid (chr22); symlink/copy.
        shutil.copyfile(bench["ref_genome"], ref_fa)
        subprocess.run([samtools, "faidx", str(ref_fa)], check=True)
    else:
        # Find the genome seqid matching the GFF seqid (exact, else first-token).
        genome_seqids = _faidx_seqids(bench["ref_genome"], samtools)
        if ref_chrom in genome_seqids:
            samtools_extract(bench["ref_genome"], [ref_chrom], ref_fa, samtools)
        else:
            raise RuntimeError(
                f"{bench['id']}: ref GFF seqid {ref_chrom!r} not found in genome "
                f"{bench['ref_genome']} (have e.g. {genome_seqids[:5]})")

    # assert harmonization
    ref_fa_seqids = set(_faidx_seqids(str(ref_fa), samtools))
    if not gff_seqids <= ref_fa_seqids:
        raise RuntimeError(f"{bench['id']}: seqid harmonization failed: "
                           f"gff={gff_seqids} not subset of fasta={ref_fa_seqids}")

    # --- target genome subset (chosen syntenic chromosome[s]) ---
    if bench.get("ready_subset"):
        shutil.copyfile(bench["tgt_genome"], tgt_fa)
        subprocess.run([samtools, "faidx", str(tgt_fa)], check=True)
        tgt_chroms = _faidx_seqids(str(tgt_fa), samtools)
    elif bench.get("tgt_chrom") == "WHOLE":
        # Very-distant pairs: genome-wide DNA synteny is gone, so there is no
        # syntenic target chromosome to pick (choose_target_chroms would minimap2
        # asm20->asm10 and then RAISE on zero hits). Feed the WHOLE target genome
        # to Liftoff + miniprot — the correct setup when the ortholog's target
        # chromosome is unknown; the reference is still subset to one chromosome
        # so the run stays tractable. Symlink (not copy) to avoid duplicating a
        # multi-GB genome; faidx the symlink so the .fai lands in the work dir.
        if tgt_fa.exists() or tgt_fa.is_symlink():
            tgt_fa.unlink()
        tgt_fa.symlink_to(os.path.abspath(bench["tgt_genome"]))
        subprocess.run([samtools, "faidx", str(tgt_fa)], check=True)
        tgt_chroms = _faidx_seqids(str(tgt_fa), samtools)
        log(f"  [synteny] WHOLE target genome: {len(tgt_chroms)} seqids "
            f"(no chromosome subsetting)")
    elif bench.get("tgt_chrom") and bench["tgt_chrom"] != "AUTO_SYNTENIC":
        # Pinned target chromosome: a concrete seqid was given (e.g. the same-species
        # giant-genome pair maize B73->Mo17, whose largest chromosome is ~313 Mb, for
        # which choose_target_chroms' whole-genome minimap2 asm20 SIGABRTs on the 2.2 Gb
        # target). Skip the synteny scan and extract the named seqid directly -- exactly
        # the chromosome AUTO_SYNTENIC would pick, for a faithful chr<->chr lift.
        # (ready_subset and WHOLE are handled above; human_mane's chr22 takes the
        # ready_subset branch, so this fires only for an explicitly pinned seqid.)
        pin = bench["tgt_chrom"]
        tgt_genome_seqids = set(_faidx_seqids(bench["tgt_genome"], samtools))
        if pin not in tgt_genome_seqids:
            raise RuntimeError(
                f"{bench['id']}: pinned tgt_chrom {pin!r} not found in target genome "
                f"{bench['tgt_genome']} (have e.g. {sorted(tgt_genome_seqids)[:5]})")
        tgt_chroms = [pin]
        samtools_extract(bench["tgt_genome"], tgt_chroms, tgt_fa, samtools)
        log(f"  [synteny] pinned target chromosome: {pin} (synteny scan skipped)")
    else:
        tgt_chroms = choose_target_chroms(ref_fa, bench["tgt_genome"],
                                          sub_dir / "synteny", minimap2, threads, log)
        tgt_genome_seqids = set(_faidx_seqids(bench["tgt_genome"], samtools))
        missing = [c for c in tgt_chroms if c not in tgt_genome_seqids]
        if missing:
            raise RuntimeError(
                f"{bench['id']}: chosen target seqid(s) {missing} not found in target "
                f"genome {bench['tgt_genome']} (have e.g. {sorted(tgt_genome_seqids)[:5]})")
        samtools_extract(bench["tgt_genome"], tgt_chroms, tgt_fa, samtools)

    # --- reference proteins for the subset ---
    n_prot = build_subset_proteins(bench, gff_info, ref_gff, ref_fa, ref_faa, log)

    manifest = {
        "id": bench["id"],
        "species": bench["species"],
        "cross_species": bench["cross_species"],
        "ref_chrom": ref_chrom,
        "tgt_chroms": tgt_chroms,
        "feature_counts": gff_info["feature_counts"],
        "n_subset_proteins": n_prot,
        "miniprot_target_space": bench.get("miniprot_target_space", "protein"),
        "protein_acc_to_mrna": gff_info["protein_acc_to_mrna"],
        "paths": {
            "ref_fa": str(ref_fa), "tgt_fa": str(tgt_fa),
            "ref_gff": str(ref_gff), "ref_faa": str(ref_faa),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    log(f"  [subset] done: {gff_info['feature_counts']}")
    return manifest


if __name__ == "__main__":  # quick manual test: python -m ... <bench_id>
    import argparse
    here = Path(__file__).resolve().parent
    reg = json.loads((here / "benchmarks.json").read_text())
    ap = argparse.ArgumentParser()
    ap.add_argument("bench_id")
    ap.add_argument("-t", "--threads", type=int, default=8)
    a = ap.parse_args()
    b = next(x for x in reg["benchmarks"] if x["id"] == a.bench_id)
    wd = here / "work" / a.bench_id
    print(json.dumps(build_subset(b, wd, reg["tools"], a.threads), indent=2)[:2000])
