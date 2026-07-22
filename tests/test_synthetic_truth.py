from __future__ import annotations

import json

from hypothesis import given, strategies as st

from benchmarks.compare.synthetic_truth import (
    build_fragmented_case,
    build_sv_case,
    deterministic_dna,
    read_fasta_sequence,
)
from benchmarks.compare.target_truth import score_target_truth
from lifton.gff3_validator import (
    validate_gff3_file,
    validate_gff3_structure,
    validate_gff3_target_bounds,
)


def _feature(gene_id, start, end, strand="+"):
    tx = gene_id + ".t"
    return (
        f"chr22\tt\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id}\n"
        f"chr22\tt\tmRNA\t{start}\t{end}\t.\t{strand}\t."
        f"\tID={tx};Parent={gene_id}\n"
        f"chr22\tt\texon\t{start}\t{end}\t.\t{strand}\t."
        f"\tID={tx}.e;Parent={tx}\n"
        f"chr22\tt\tCDS\t{start}\t{end}\t.\t{strand}\t0"
        f"\tID={tx}.c;Parent={tx}\n"
    )


def _inputs(tmp_path, length, body):
    fasta = tmp_path / "source.fa"
    gff = tmp_path / "source.gff3"
    sequence = ("ACGTTGCA" * ((length + 7) // 8))[:length]
    fasta.write_text(">chr22 source\n" + sequence + "\n")
    gff.write_text("##gff-version 3\n" + body)
    return fasta, gff, sequence


def test_fragmented_case_has_exact_coordinates_and_golden_output(tmp_path):
    fasta, gff, sequence = _inputs(
        tmp_path,
        80,
        _feature("g1", 5, 15)
        + _feature("crossing", 18, 25)
        + _feature("g2", 55, 65, "-"),
    )

    outputs = build_fragmented_case(
        fasta, gff, tmp_path / "fragmented", cuts=(20, 50),
    )

    target = outputs["target_fasta"]
    assert read_fasta_sequence(target, "chr22_frag_001") == sequence[:20]
    assert read_fasta_sequence(target, "chr22_frag_002") == sequence[20:50]
    assert read_fasta_sequence(target, "chr22_frag_003") == sequence[50:]
    truth_text = outputs["truth_gff"].read_text()
    assert (
        "chr22_frag_001\tLiftOn-synthetic-truth\tgene\t5\t15"
        in truth_text
    )
    assert (
        "chr22_frag_003\tLiftOn-synthetic-truth\tgene\t5\t15\t.\t-"
        in truth_text
    )
    assert "ID=crossing" not in truth_text
    mapping = json.loads(outputs["ortholog_map"].read_text())
    crossing = next(
        row for row in mapping["mappings"]
        if row["source_id"] == "crossing"
    )
    assert crossing == {
        "feature_type": "gene",
        "source_id": "crossing",
        "status": "breakpoint_crossing",
        "truth_ids": [],
    }
    assert validate_gff3_target_bounds(
        str(outputs["truth_gff"]), str(target),
    ).is_valid
    assert validate_gff3_structure(str(outputs["truth_gff"])).is_valid
    assert validate_gff3_file(str(outputs["truth_gff"])).is_valid


def test_sv_case_maps_delete_invert_duplicate_insert_exactly(tmp_path):
    fasta, gff, _sequence = _inputs(
        tmp_path,
        100,
        _feature("before", 2, 5)
        + _feature("deleted", 12, 18)
        + _feature("deletion_crossing", 8, 13)
        + _feature("inverted", 32, 38)
        + _feature("duplicated", 52, 58)
        + _feature("insertion_crossing", 68, 74)
        + _feature("after", 75, 80),
    )

    outputs = build_sv_case(
        fasta,
        gff,
        tmp_path / "sv",
        deletion=(11, 20),
        inversion=(31, 40),
        duplication=(51, 60),
        insert_after=70,
        insertion_length=8,
        seed="golden",
    )

    assert len(read_fasta_sequence(outputs["target_fasta"], "chr22")) == 108
    truth = outputs["truth_gff"].read_text()
    assert (
        "chr22\tLiftOn-synthetic-truth\tgene\t23\t29\t.\t-"
        "\t.\tID=inverted"
    ) in truth
    assert (
        "chr22\tLiftOn-synthetic-truth\tgene\t42\t48\t.\t+\t."
        "\tID=duplicated__synthetic_copy1"
    ) in truth
    assert (
        "chr22\tLiftOn-synthetic-truth\tgene\t52\t58\t.\t+\t."
        "\tID=duplicated__synthetic_copy2"
    ) in truth
    assert (
        "chr22\tLiftOn-synthetic-truth\tgene\t83\t88\t.\t+\t.\tID=after"
    ) in truth
    assert "ID=deleted" not in truth
    assert "ID=deletion_crossing" not in truth
    assert "ID=insertion_crossing" not in truth
    assert validate_gff3_target_bounds(
        str(outputs["truth_gff"]), str(outputs["target_fasta"]),
    ).is_valid
    assert validate_gff3_structure(str(outputs["truth_gff"])).is_valid
    assert validate_gff3_file(str(outputs["truth_gff"])).is_valid

    mapping = json.loads(outputs["ortholog_map"].read_text())
    mappings = {row["source_id"]: row for row in mapping["mappings"]}
    assert mappings["deleted"]["status"] == "deleted"
    assert mappings["deletion_crossing"]["status"] == "breakpoint_crossing"
    assert mappings["duplicated"]["truth_ids"] == [
        "duplicated__synthetic_copy1",
        "duplicated__synthetic_copy2",
    ]

    prediction = tmp_path / "synthetic-perfect-prediction.gff3"
    prediction.write_text(
        truth.replace("__synthetic_copy1", "").replace(
            "__synthetic_copy2", "_1",
        )
    )
    metrics = score_target_truth(
        prediction,
        outputs["truth_gff"],
        ortholog_map=outputs["ortholog_map"],
    )
    assert metrics["gene"]["locus"]["f1"] == 1.0
    assert metrics["gene"]["copy"]["f1"] == 1.0
    assert metrics["transcript"]["strand"]["f1"] == 1.0
    assert metrics["structure"]["intron_chain"]["f1"] == 1.0
    assert metrics["structure"]["exon"]["f1"] == 1.0
    assert metrics["structure"]["CDS"]["f1"] == 1.0


def test_synthetic_manifests_are_reproducible_across_output_directories(
        tmp_path):
    fasta, gff, _sequence = _inputs(
        tmp_path,
        80,
        _feature("g1", 5, 15) + _feature("g2", 55, 65),
    )

    first = build_fragmented_case(
        fasta, gff, tmp_path / "first", cuts=(20, 50),
    )
    second = build_fragmented_case(
        fasta, gff, tmp_path / "second", cuts=(20, 50),
    )

    for key in ("target_fasta", "truth_gff", "ortholog_map", "manifest"):
        assert first[key].read_bytes() == second[key].read_bytes()


def test_synthetic_mapping_uses_decoded_logical_gff3_ids(tmp_path):
    body = (
        "chr22\tt\tgene\t5\t15\t.\t+\t.\tID=gene%3Ag%2Cone\n"
        "chr22\tt\tmRNA\t5\t15\t.\t+\t."
        "\tID=tx%3Aone;Parent=gene%3Ag%2Cone\n"
        "chr22\tt\texon\t5\t15\t.\t+\t."
        "\tID=exon%3Aone;Parent=tx%3Aone\n"
        "chr22\tt\tCDS\t5\t15\t.\t+\t0"
        "\tID=cds%3Aone;Parent=tx%3Aone\n"
    )
    fasta, gff, _sequence = _inputs(tmp_path, 40, body)

    outputs = build_fragmented_case(
        fasta, gff, tmp_path / "encoded", cuts=(20,),
    )

    mapping = json.loads(outputs["ortholog_map"].read_text())
    by_source = {row["source_id"]: row for row in mapping["mappings"]}
    assert by_source["gene:g,one"]["truth_ids"] == ["gene:g,one"]
    assert by_source["tx:one"]["truth_ids"] == ["tx:one"]
    metrics = score_target_truth(
        outputs["truth_gff"],
        outputs["truth_gff"],
        ortholog_map=outputs["ortholog_map"],
    )
    assert metrics["gene"]["locus"]["f1"] == 1.0
    assert metrics["transcript"]["locus"]["f1"] == 1.0


@given(
    length=st.integers(min_value=0, max_value=500),
    seed=st.text(min_size=0, max_size=30),
)
def test_deterministic_dna_property(length, seed):
    first = deterministic_dna(length, seed)
    second = deterministic_dna(length, seed)

    assert first == second
    assert len(first) == length
    assert set(first) <= set("ACGT")
