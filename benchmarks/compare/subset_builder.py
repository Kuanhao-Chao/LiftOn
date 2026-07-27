"""Subset builder: reduce each benchmark to ONE chromosome pair.

Produces, under ``work/<id>/subset/``:
  ref.chrom.fa / .fai      subset reference genome (seqid harmonized to GFF)
  tgt.chrom.fa / .fai      subset target genome (chosen syntenic chromosome[s])
  ref.chrom.gff3           subset reference annotation (single ref chromosome)
  ref.proteins.subset.faa  reference proteins restricted to the subset
  subset.manifest.json     chosen chroms, header maps, counts
  synteny/ref_chrom_to_tgt.paf, target_chroms.txt   (cross-species / AUTO)

Reuse is fail-closed: a versioned manifest binds source/tool/builder hashes,
parameters, commands, and every emitted artifact.  The legacy ``.done`` file is
written for compatibility but is never trusted as cache evidence.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence, TextIO


SUBSET_MANIFEST_SCHEMA_VERSION = 2
SUBSET_BUILDER_VERSION = "content-addressed-subset-v2"
SUBSET_MANIFEST_KEYS = {
    "schema_version",
    "builder_version",
    "request_sha256",
    "request",
    "id",
    "species",
    "cross_species",
    "ref_chrom",
    "tgt_chroms",
    "feature_counts",
    "n_subset_proteins",
    "miniprot_target_space",
    "annotation_format",
    "protein_acc_to_mrna",
    "paths",
    "commands",
    "outputs",
}
_GTF_ATTRIBUTE_RE = re.compile(
    r'''^\s*([^\s;]+)\s+(?:"([^"]*)"|([^;\s]+))\s*$'''
)


# ---------------------------------------------------------------------------
# small GFF streaming helpers (no gffutils DB build needed for subsetting)
# ---------------------------------------------------------------------------

def _open_text(path: str | Path) -> TextIO:
    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="strict")


def _gff3_attrs(col9: str) -> dict[str, str]:
    """Parse a GFF3 column-9 attribute string into a dict (first value wins)."""
    result = {}
    for kv in col9.rstrip(";").split(";"):
        kv = kv.strip()
        if not kv or "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        result.setdefault(k, v)
    return result


def _gtf_attrs(col9: str) -> dict[str, str]:
    """Parse GTF attributes without silently treating them as GFF3."""

    result = {}
    for raw in col9.split(";"):
        raw = raw.strip()
        if not raw:
            continue
        match = _GTF_ATTRIBUTE_RE.fullmatch(raw)
        if match is None:
            raise RuntimeError(f"malformed GTF attribute: {raw!r}")
        result.setdefault(
            match.group(1),
            match.group(2) if match.group(2) is not None else match.group(3),
        )
    return result


def _attrs(col9: str, annotation_format: str = "GFF3") -> dict[str, str]:
    """Parse column 9 according to the declared annotation dialect."""

    return (
        _gtf_attrs(col9)
        if annotation_format.upper() == "GTF"
        else _gff3_attrs(col9)
    )


def detect_annotation_format(path: str | Path, declared: str | None = None) -> str:
    """Return ``GTF`` or ``GFF3`` and reject a conflicting declaration."""

    declared_format = str(declared or "").upper() or None
    if declared_format not in {None, "GFF3", "GTF"}:
        raise RuntimeError(f"unsupported annotation format: {declared!r}")
    observed = None
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise RuntimeError(
                    f"{path}: line {line_number} has {len(columns)} columns"
                )
            attributes = columns[8].strip()
            observed = (
                "GFF3"
                if attributes in {"", "."}
                or re.match(r"^[^;\s=]+\s*=", attributes)
                else "GTF"
            )
            break
    if observed is None:
        raise RuntimeError(f"annotation contains no feature records: {path}")
    if declared_format is not None and declared_format != observed:
        raise RuntimeError(
            f"declared annotation format {declared_format} differs from "
            f"observed {observed}: {path}"
        )
    return observed


def gff_mrna_counts_per_seqid(gff_path: str) -> Counter:
    """Count coding transcript rows per seqid in GFF3 or GTF."""
    counts: Counter = Counter()
    with _open_text(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            if parts[2].lower() in {"mrna", "transcript"}:
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


def subset_gff(
    gff_path: str,
    keep_seqid: str,
    out_path: Path,
    *,
    annotation_format: str | None = None,
) -> dict:
    """Subset GFF3 or GTF while preserving the source annotation dialect."""

    annotation_format = detect_annotation_format(gff_path, annotation_format)
    counts: Counter = Counter()
    protein_accs: set[str] = set()
    acc_to_mrna: dict[str, str] = {}
    cds_parent_protein: list[tuple[str, str]] = []  # (parent_mrna, protein_acc)
    with _open_text(gff_path) as fh, open(
        out_path, "w", encoding="utf-8", newline="\n",
    ) as out:
        if annotation_format == "GFF3":
            out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#"):
                if annotation_format == "GTF":
                    out.write(line if line.endswith("\n") else line + "\n")
                elif line.startswith("##sequence-region"):
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
                a = _attrs(parts[8], annotation_format)
                pacc = a.get("protein_id") or a.get("Name")
                parent = (
                    a.get("transcript_id", "")
                    if annotation_format == "GTF"
                    else a.get("Parent", "").split(",", 1)[0]
                )
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
        "annotation_format": annotation_format,
    }


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(document: Any) -> str:
    payload = json.dumps(
        document,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _file_record(path: str | Path) -> dict[str, Any]:
    absolute = Path(os.path.abspath(Path(path).expanduser()))
    if not absolute.is_file() or absolute.stat().st_size <= 0:
        raise RuntimeError(f"cache artifact is missing or empty: {absolute}")
    return {
        "path": str(absolute),
        "size": absolute.stat().st_size,
        "sha256": _sha256_file(absolute),
    }


def _resolve_tool(executable: str | Path) -> Path:
    text = str(executable)
    candidate = (
        shutil.which(text)
        if os.sep not in text and (os.altsep is None or os.altsep not in text)
        else text
    )
    if candidate is None:
        raise RuntimeError(f"required executable is not on PATH: {text}")
    path = Path(candidate).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise RuntimeError(f"required executable is not executable: {path}")
    return path


def _tool_record(executable: str | Path) -> dict[str, Any]:
    resolved = _resolve_tool(executable)
    return {
        "configured": str(executable),
        "executable": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": _sha256_file(resolved),
    }


def _subset_request(
    bench: Mapping[str, Any],
    tools: Mapping[str, Any],
    *,
    threads: int,
    annotation_format: str,
) -> dict[str, Any]:
    target_chromosome = bench.get("tgt_chrom") or "AUTO_SYNTENIC"
    inputs = {
        role: _file_record(bench[role])
        for role in ("ref_genome", "ref_gff", "tgt_genome")
    }
    if bench.get("ref_proteins"):
        inputs["ref_proteins"] = _file_record(bench["ref_proteins"])
    tool_records = {
        "samtools": _tool_record(tools["samtools_bin"]),
    }
    if (
        not bench.get("ready_subset")
        and target_chromosome == "AUTO_SYNTENIC"
    ):
        tool_records["minimap2"] = _tool_record(tools["minimap2_bin"])
    request = {
        "schema_version": SUBSET_MANIFEST_SCHEMA_VERSION,
        "builder": {
            "id": SUBSET_BUILDER_VERSION,
            "source_sha256": _sha256_file(Path(__file__).resolve()),
        },
        "benchmark": {
            "id": bench["id"],
            "species": bench["species"],
            "cross_species": bool(bench["cross_species"]),
            "ref_chrom": bench["ref_chrom"],
            "tgt_chrom": target_chromosome,
            "ready_subset": bool(bench.get("ready_subset")),
            "miniprot_target_space": bench.get(
                "miniprot_target_space", "protein"
            ),
            "annotation_format": annotation_format,
            "threads": threads,
        },
        "inputs": inputs,
        "tools": tool_records,
    }
    request["request_sha256"] = _canonical_sha256(request)
    return request


def _expected_cached_paths(
    manifest_path: Path,
    request: Mapping[str, Any],
) -> tuple[dict[str, str], dict[str, str]]:
    """Return the only logical artifact paths valid for this cache request."""

    sub_dir = manifest_path.parent
    extension = (
        "gtf"
        if request["benchmark"]["annotation_format"] == "GTF"
        else "gff3"
    )
    paths = {
        "ref_fa": str(Path(os.path.abspath(sub_dir / "ref.chrom.fa"))),
        "tgt_fa": str(Path(os.path.abspath(sub_dir / "tgt.chrom.fa"))),
        "ref_gff": str(Path(os.path.abspath(
            sub_dir / f"ref.chrom.{extension}"
        ))),
        "ref_faa": str(Path(os.path.abspath(
            sub_dir / "ref.proteins.subset.faa"
        ))),
    }
    outputs = {
        **paths,
        "ref_fai": str(Path(os.path.abspath(paths["ref_fa"] + ".fai"))),
        "tgt_fai": str(Path(os.path.abspath(paths["tgt_fa"] + ".fai"))),
    }
    benchmark = request["benchmark"]
    if (
        not benchmark["ready_subset"]
        and benchmark["tgt_chrom"] == "AUTO_SYNTENIC"
    ):
        synteny_dir = sub_dir / "synteny"
        outputs.update({
            "synteny_paf": str(Path(os.path.abspath(
                synteny_dir / "ref_chrom_to_tgt.paf"
            ))),
            "target_chroms": str(Path(os.path.abspath(
                synteny_dir / "target_chroms.txt"
            ))),
        })
    return paths, outputs


def _cached_subset_manifest(
    manifest_path: Path,
    request: Mapping[str, Any],
) -> dict[str, Any] | None:
    try:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            not isinstance(document, dict)
            or set(document) != SUBSET_MANIFEST_KEYS
            or document.get("schema_version") != SUBSET_MANIFEST_SCHEMA_VERSION
            or document.get("builder_version") != SUBSET_BUILDER_VERSION
            or document.get("request") != request
            or document.get("request_sha256") != request["request_sha256"]
            or not isinstance(document.get("outputs"), dict)
        ):
            return None
        benchmark = request["benchmark"]
        if (
            document["id"] != benchmark["id"]
            or document["species"] != benchmark["species"]
            or document["cross_species"] != benchmark["cross_species"]
            or document["miniprot_target_space"]
            != benchmark["miniprot_target_space"]
            or document["annotation_format"] != benchmark["annotation_format"]
        ):
            return None
        expected_paths, expected_outputs = _expected_cached_paths(
            manifest_path, request,
        )
        if (
            document["paths"] != expected_paths
            or set(document["outputs"]) != set(expected_outputs)
        ):
            return None
        for label, expected_path in expected_outputs.items():
            record = document["outputs"][label]
            if (
                not isinstance(record, dict)
                or set(record) != {"path", "size", "sha256"}
                or record["path"] != expected_path
                or _file_record(record["path"]) != record
            ):
                return None
        if not isinstance(document["commands"], list) or any(
            not isinstance(command, dict)
            or command.get("shell") is not False
            or not isinstance(command.get("argv"), list)
            or not command["argv"]
            for command in document["commands"]
        ):
            return None
        return document
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def verify_cached_subset(
    bench: Mapping[str, Any],
    work_dir: str | Path,
    tools: Mapping[str, Any],
    *,
    threads: int,
) -> bool:
    """Verify a subset cache against its current inputs, tools, and builder."""

    try:
        annotation_format = detect_annotation_format(
            bench["ref_gff"], bench.get("input_format"),
        )
        request = _subset_request(
            bench,
            tools,
            threads=threads,
            annotation_format=annotation_format,
        )
        manifest_path = Path(work_dir) / "subset" / "subset.manifest.json"
        return _cached_subset_manifest(manifest_path, request) is not None
    except (KeyError, OSError, RuntimeError, TypeError, ValueError):
        return False


def _atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


# ---------------------------------------------------------------------------
# FASTA helpers (samtools faidx)
# ---------------------------------------------------------------------------

def _command_record(
    commands: list[dict[str, Any]] | None,
    label: str,
    argv: Sequence[str | Path],
    *,
    stdout: str | Path | None = None,
) -> None:
    if commands is None:
        return
    record = {
        "label": label,
        "argv": [str(value) for value in argv],
        "shell": False,
    }
    if stdout is not None:
        record["stdout"] = str(Path(stdout).resolve())
    commands.append(record)


def _faidx_seqids(
    fa: str,
    samtools: str,
    commands: list[dict[str, Any]] | None = None,
) -> list[str]:
    fai = fa + ".fai"
    if not os.path.exists(fai):
        argv = [samtools, "faidx", fa]
        _command_record(commands, "index_fasta", argv)
        subprocess.run(argv, check=True)
    out = []
    with open(fai) as fh:
        for line in fh:
            out.append(line.split("\t")[0])
    return out


def samtools_extract(
    genome: str,
    seqids: list[str],
    out_fa: Path,
    samtools: str,
    rename_to: str | None = None,
    commands: list[dict[str, Any]] | None = None,
) -> None:
    """Extract ``seqids`` from ``genome`` into ``out_fa`` and index it.
    If ``rename_to`` is given (single seqid), the output header is rewritten."""
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    extract_argv = [samtools, "faidx", genome, *seqids]
    _command_record(
        commands, "extract_fasta_sequences", extract_argv, stdout=out_fa,
    )
    with open(out_fa, "w") as out:
        proc = subprocess.run(extract_argv,
                              stdout=subprocess.PIPE, check=True, text=True)
        text = proc.stdout
        if rename_to is not None and len(seqids) == 1:
            lines = text.splitlines(keepends=True)
            if lines and lines[0].startswith(">"):
                lines[0] = f">{rename_to}\n"
            text = "".join(lines)
        out.write(text)
    output_index = Path(str(out_fa) + ".fai")
    if output_index.exists():
        output_index.unlink()
    index_argv = [samtools, "faidx", str(out_fa)]
    _command_record(commands, "index_subset_fasta", index_argv)
    subprocess.run(index_argv, check=True)


# ---------------------------------------------------------------------------
# synteny: pick target chromosome(s) for the chosen reference chromosome
# ---------------------------------------------------------------------------

def choose_target_chroms(ref_chrom_fa: Path, tgt_genome: str, out_dir: Path,
                         minimap2: str, threads: int, log=print,
                         commands: list[dict[str, Any]] | None = None) -> list[str]:
    """minimap2-align the subset reference chromosome to the full target genome,
    return the dominant syntenic target seqid(s) (>= max(10kb, 10% of matches))."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paf = out_dir / "ref_chrom_to_tgt.paf"

    def _run_minimap(preset: str) -> dict:
        argv = [minimap2, "-x", preset, "-t", str(threads),
                "--secondary=no", tgt_genome, str(ref_chrom_fa)]
        _command_record(
            commands, f"target_synteny_{preset}", argv, stdout=paf,
        )
        with open(paf, "w") as fh:
            subprocess.run(argv,
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
                    key=lambda s: (-matches[s], s))
    if not chosen:  # fall back to the single best
        chosen = [min(matches, key=lambda value: (-matches[value], value))]
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
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise RuntimeError("subset threads must be a positive integer")

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

    annotation_format = detect_annotation_format(
        bench["ref_gff"], bench.get("input_format"),
    )
    request = _subset_request(
        bench,
        tools,
        threads=threads,
        annotation_format=annotation_format,
    )
    if not force and manifest_path.is_file():
        cached = _cached_subset_manifest(manifest_path, request)
        if cached is not None:
            log(f"  [subset] verified cache ({bench['id']})")
            return cached

    sub_dir.mkdir(parents=True, exist_ok=True)
    samtools = request["tools"]["samtools"]["executable"]
    minimap2_record = request["tools"].get("minimap2")
    minimap2 = (
        minimap2_record["executable"]
        if minimap2_record is not None
        else str(tools["minimap2_bin"])
    )
    commands: list[dict[str, Any]] = []

    ref_chrom = choose_ref_chrom(bench)
    log(f"  [subset] ref chromosome = {ref_chrom}")

    ref_fa = sub_dir / "ref.chrom.fa"
    tgt_fa = sub_dir / "tgt.chrom.fa"
    ref_gff = sub_dir / (
        "ref.chrom.gtf" if annotation_format == "GTF" else "ref.chrom.gff3"
    )
    ref_faa = sub_dir / "ref.proteins.subset.faa"

    # --- subset reference annotation (single chromosome, dialect preserved) ---
    gff_info = subset_gff(
        bench["ref_gff"],
        ref_chrom,
        ref_gff,
        annotation_format=annotation_format,
    )
    gff_seqids = {ref_chrom}

    # --- subset reference genome FASTA, harmonize header to the GFF seqid ---
    if bench.get("ready_subset"):
        # Ready chr22 FASTAs already use the GFF seqid (chr22); symlink/copy.
        shutil.copyfile(bench["ref_genome"], ref_fa)
        ref_index = Path(str(ref_fa) + ".fai")
        if ref_index.exists():
            ref_index.unlink()
        argv = [samtools, "faidx", str(ref_fa)]
        _command_record(commands, "index_ready_reference_fasta", argv)
        subprocess.run(argv, check=True)
    else:
        # Find the genome seqid matching the GFF seqid (exact, else first-token).
        genome_seqids = _faidx_seqids(
            bench["ref_genome"], samtools, commands,
        )
        if ref_chrom in genome_seqids:
            samtools_extract(
                bench["ref_genome"], [ref_chrom], ref_fa, samtools,
                commands=commands,
            )
        else:
            raise RuntimeError(
                f"{bench['id']}: ref GFF seqid {ref_chrom!r} not found in genome "
                f"{bench['ref_genome']} (have e.g. {genome_seqids[:5]})")

    # assert harmonization
    ref_fa_seqids = set(_faidx_seqids(str(ref_fa), samtools, commands))
    if not gff_seqids <= ref_fa_seqids:
        raise RuntimeError(f"{bench['id']}: seqid harmonization failed: "
                           f"gff={gff_seqids} not subset of fasta={ref_fa_seqids}")

    # --- target genome subset (chosen syntenic chromosome[s]) ---
    if bench.get("ready_subset"):
        shutil.copyfile(bench["tgt_genome"], tgt_fa)
        target_index = Path(str(tgt_fa) + ".fai")
        if target_index.exists():
            target_index.unlink()
        argv = [samtools, "faidx", str(tgt_fa)]
        _command_record(commands, "index_ready_target_fasta", argv)
        subprocess.run(argv, check=True)
        tgt_chroms = _faidx_seqids(str(tgt_fa), samtools, commands)
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
        target_index = Path(str(tgt_fa) + ".fai")
        if target_index.exists():
            target_index.unlink()
        argv = [samtools, "faidx", str(tgt_fa)]
        _command_record(commands, "index_whole_target_fasta", argv)
        subprocess.run(argv, check=True)
        tgt_chroms = _faidx_seqids(str(tgt_fa), samtools, commands)
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
        tgt_genome_seqids = set(
            _faidx_seqids(bench["tgt_genome"], samtools, commands)
        )
        if pin not in tgt_genome_seqids:
            raise RuntimeError(
                f"{bench['id']}: pinned tgt_chrom {pin!r} not found in target genome "
                f"{bench['tgt_genome']} (have e.g. {sorted(tgt_genome_seqids)[:5]})")
        tgt_chroms = [pin]
        samtools_extract(
            bench["tgt_genome"], tgt_chroms, tgt_fa, samtools,
            commands=commands,
        )
        log(f"  [synteny] pinned target chromosome: {pin} (synteny scan skipped)")
    else:
        tgt_chroms = choose_target_chroms(ref_fa, bench["tgt_genome"],
                                          sub_dir / "synteny", minimap2, threads, log,
                                          commands)
        tgt_genome_seqids = set(
            _faidx_seqids(bench["tgt_genome"], samtools, commands)
        )
        missing = [c for c in tgt_chroms if c not in tgt_genome_seqids]
        if missing:
            raise RuntimeError(
                f"{bench['id']}: chosen target seqid(s) {missing} not found in target "
                f"genome {bench['tgt_genome']} (have e.g. {sorted(tgt_genome_seqids)[:5]})")
        samtools_extract(
            bench["tgt_genome"], tgt_chroms, tgt_fa, samtools,
            commands=commands,
        )

    # --- reference proteins for the subset ---
    n_prot = build_subset_proteins(bench, gff_info, ref_gff, ref_fa, ref_faa, log)

    paths, output_paths = _expected_cached_paths(manifest_path, request)
    manifest = {
        "schema_version": SUBSET_MANIFEST_SCHEMA_VERSION,
        "builder_version": SUBSET_BUILDER_VERSION,
        "request_sha256": request["request_sha256"],
        "request": request,
        "id": bench["id"],
        "species": bench["species"],
        "cross_species": bench["cross_species"],
        "ref_chrom": ref_chrom,
        "tgt_chroms": tgt_chroms,
        "feature_counts": gff_info["feature_counts"],
        "n_subset_proteins": n_prot,
        "miniprot_target_space": bench.get("miniprot_target_space", "protein"),
        "annotation_format": annotation_format,
        "protein_acc_to_mrna": gff_info["protein_acc_to_mrna"],
        "paths": paths,
        "commands": commands,
        "outputs": {
            label: _file_record(path)
            for label, path in sorted(output_paths.items())
        },
    }
    _atomic_write_json(manifest_path, manifest)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"{request['request_sha256']}\n", encoding="ascii")
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
