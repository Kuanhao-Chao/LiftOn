"""Pin deterministic, genomic-order extra-copy numbering (FIX 1).

The vendored Liftoff numbered gene copies in alignment-encounter order
(`write_new_gff.write_new_gff` sorted parents by `.id` only, a stable sort that
preserved the non-deterministic dict-insertion order), so the `_N` copy suffix on
a given genomic locus could differ run-to-run. The fix sorts parents by
`(seqid, start, end, id)` before `finalize_parent_features` assigns copy_num, so
`_0` is always the 5'-most copy on each sequence and the numbering is
deterministic regardless of the order copies were discovered.

These tests exercise the REAL `write_new_gff` end-to-end (so removing the sort
fails the test), using gffutils.Feature objects and parsing the emitted GFF3.
"""
import types

import gffutils
import pytest

from lifton.liftoff import write_new_gff


def _feat(ftype, seqid, start, end, fid, parent=None, copy_id=None):
    attrs = {"ID": [fid]}
    if parent is not None:
        attrs["Parent"] = [parent]
    if ftype == "gene":
        attrs["copy_id"] = [copy_id or fid]
        attrs["coverage"] = ["1.0"]
        attrs["sequence_ID"] = ["1.0"]
    return gffutils.Feature(
        seqid=seqid, source="Liftoff", featuretype=ftype, id=fid,
        start=start, end=end, score=".", strand="+", frame=".", attributes=attrs)


def _one_copy(seqid, start):
    """A gene→mRNA→exon→CDS stack for one copy, keyed by a unique copy_id."""
    copy_id = f"gene1@{seqid}:{start}"
    end = start + 500
    gene = _feat("gene", seqid, start, end, "gene1", copy_id=copy_id)
    mrna = _feat("mRNA", seqid, start, end, "rna1", parent="gene1")
    exon = _feat("exon", seqid, start, end, "exon1", parent="rna1")
    cds = _feat("CDS", seqid, start, end, "cds1", parent="rna1")
    return copy_id, [gene, mrna, exon, cds]


def _run_write(lifted_features, tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    out = tmp_path / "liftoff.gff3"
    args = types.SimpleNamespace(output=str(out), a=0.5, s=0.5)
    feature_db = types.SimpleNamespace(dialect={"fmt": "gff3"})
    write_new_gff.write_new_gff(lifted_features, args, feature_db)
    # parse: map gene-line start -> emitted gene ID
    start_to_id = {}
    for line in out.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if f[2] == "gene":
            gid = next(a[3:] for a in f[8].split(";") if a.startswith("ID="))
            start_to_id[int(f[3])] = gid
    return start_to_id


def test_copy_numbering_is_genomic_order(tmp_path):
    """copy_0 (bare id) is the 5'-most copy; _1, _2 ascend by coordinate."""
    # dict insertion order is NON-genomic (5000, 1000, 3000) — as if discovered
    # by the aligner out of order.
    lifted = {}
    for s in (5000, 1000, 3000):
        cid, feats = _one_copy("chr1", s)
        lifted[cid] = feats
    start_to_id = _run_write(lifted, tmp_path)
    assert start_to_id[1000] == "gene1"      # 5'-most → copy 0 (no suffix)
    assert start_to_id[3000] == "gene1_1"
    assert start_to_id[5000] == "gene1_2"


def test_copy_numbering_is_deterministic_across_input_order(tmp_path):
    """A different discovery order yields the SAME genomic numbering."""
    order_a, order_b = (5000, 1000, 3000), (3000, 5000, 1000)
    res = []
    for order in (order_a, order_b):
        lifted = {}
        for s in order:
            cid, feats = _one_copy("chr1", s)
            lifted[cid] = feats
        res.append(_run_write(lifted, tmp_path / f"o{order[0]}"))
    assert res[0] == res[1]
    assert res[0] == {1000: "gene1", 3000: "gene1_1", 5000: "gene1_2"}


def test_copy_numbering_per_sequence(tmp_path):
    """Numbering restarts per (gene); copies on different seqids order by (seqid,start)."""
    lifted = {}
    for seqid, s in [("chr2", 200), ("chr1", 900), ("chr1", 100)]:
        cid, feats = _one_copy(seqid, s)
        lifted[cid] = feats
    start_to_id = _run_write(lifted, tmp_path)
    # global genomic order: chr1:100, chr1:900, chr2:200  → _0, _1, _2
    assert start_to_id[100] == "gene1"
    assert start_to_id[900] == "gene1_1"
    assert start_to_id[200] == "gene1_2"


def test_child_ids_carry_the_copy_suffix(tmp_path):
    """A non-zero copy's child features (mRNA/exon/CDS) get the `_N` suffix."""
    lifted = {}
    for s in (5000, 1000):
        cid, feats = _one_copy("chr1", s)
        lifted[cid] = feats
    out = tmp_path / "liftoff.gff3"
    args = types.SimpleNamespace(output=str(out), a=0.5, s=0.5)
    feature_db = types.SimpleNamespace(dialect={"fmt": "gff3"})
    write_new_gff.write_new_gff(lifted, args, feature_db)
    text = out.read_text()
    # the 5000 copy is copy_1 → its mRNA must be rna1_1 with Parent gene1_1
    assert "ID=rna1_1;Parent=gene1_1" in text
    # the 1000 copy is copy_0 → bare ids
    assert "ID=rna1;Parent=gene1" in text


def test_cross_sequence_child_branch_uses_matching_parent_copy(tmp_path):
    """A trans-spliced branch must not attach to a different-seqid copy."""
    first_key = "gene1@chr1:100"
    second_key = "gene1@chr2:100"
    first_parent = _feat("gene", "chr1", 100, 200, "gene1", copy_id=first_key)
    second_parent = _feat("gene", "chr2", 100, 200, "gene1", copy_id=second_key)
    transcript = _feat("mRNA", "chr2", 100, 200, "rna1", parent="gene1")
    exon = _feat("exon", "chr2", 100, 200, "exon1", parent="rna1")
    lifted = {
        first_key: [first_parent],
        second_key: [second_parent],
    }
    lifted[first_key].extend([transcript, exon])

    out = tmp_path / "trans-spliced.gff3"
    write_new_gff.write_new_gff(
        lifted,
        types.SimpleNamespace(output=str(out), a=0.5, s=0.5),
        types.SimpleNamespace(dialect={"fmt": "gff3"}),
    )

    lines = [line for line in out.read_text().splitlines() if line and not line.startswith("#")]
    transcript_line = next(line for line in lines if "\tmRNA\t" in line)
    assert transcript_line.split("\t")[0] == "chr2"
    assert "ID=rna1_1;Parent=gene1_1" in transcript_line
    assert not any("ID=rna1;Parent=gene1" in line for line in lines)
