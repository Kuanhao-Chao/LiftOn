from __future__ import annotations

import gzip

import pytest

from benchmarks.compare import provider_orthologs


def test_gene_id_index_excludes_ambiguous_models(tmp_path):
    annotation = tmp_path / "models.gff3"
    annotation.write_text(
        "##gff-version 3\n"
        "chr1\tNCBI\tgene\t1\t10\t.\t+\t.\tID=g1;Dbxref=GeneID:10\n"
        "chr1\tNCBI\tgene\t20\t30\t.\t+\t.\tID=g2;Dbxref=GeneID:20\n"
        "chr1\tNCBI\tgene\t40\t50\t.\t+\t.\tID=g3;Dbxref=GeneID:20\n"
        "chr1\tNCBI\tgene\t60\t70\t.\t+\t.\tID=g4;Dbxref=GeneID:30,GeneID:31\n",
        encoding="utf-8",
    )

    assert provider_orthologs.gene_id_index(annotation) == {"10": "g1"}


def test_relationship_parser_orients_requested_taxon_pair(tmp_path):
    table = tmp_path / "orthologs.gz"
    with gzip.open(table, "wt", encoding="utf-8") as handle:
        handle.write(
            "#tax_id\tGeneID\trelationship\tOther_tax_id\tOther_GeneID\n"
            "9606\t1\tOrtholog\t9595\t11\n"
            "9595\t12\tOrtholog\t9606\t2\n"
            "9606\t3\tOther\t9595\t13\n"
        )

    parsed = provider_orthologs.relevant_relationships(
        table, {"human_gorilla": (9606, 9595)},
    )

    assert dict(parsed["human_gorilla"]) == {
        "1": {"11"}, "2": {"12"},
    }


def test_provider_mapping_keeps_only_one_to_one_annotation_groups():
    mapping = provider_orthologs.build_mapping(
        "pair",
        {"1": {"11"}, "2": {"12"}, "3": {"12"}, "4": {"14", "15"}},
        {"1": "sg1", "2": "sg2", "3": "sg3", "4": "sg4"},
        {"11": "tg11", "12": "tg12", "14": "tg14", "15": "tg15"},
    )

    assert mapping["metadata"]["groups"] == 1
    assert mapping["mappings"] == [{
        "source_id": "sg1",
        "truth_ids": ["tg11"],
        "feature_type": "gene",
        "status": "retained",
        "evidence": {
            "source_gene_id": "1",
            "target_gene_id": "11",
            "relationship": "Ortholog",
        },
    }]


def test_relationship_parser_rejects_wrong_header(tmp_path):
    table = tmp_path / "bad.gz"
    with gzip.open(table, "wt", encoding="utf-8") as handle:
        handle.write("#wrong\theader\n")

    with pytest.raises(provider_orthologs.ProviderOrthologError, match="header"):
        provider_orthologs.relevant_relationships(
            table, {"pair": (1, 2)},
        )
