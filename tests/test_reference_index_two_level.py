"""A gene that has direct exons AND a transcript must still index the transcript.

`lifton_utils.get_ref_liffover_features` branches on whether a top-level
feature has direct level-1 exon children:

    if len(exon_children) > 0:
        __process_ref_liffover_features(locus, ref_db, None)   # indexes nothing
    else:
        for transcript in children(locus, level=1):
            ...  ref_features_reverse_dict[transcript.id] = locus.id  ...

`ref_features_reverse_dict` and `ref_trans_exon_num_dict` are therefore built
only for 3-level `gene -> transcript -> exon` annotations. RefSeq's organellar
convention writes a plastid gene's exons **twice** -- once directly under the
gene and once under its mRNA -- so such a gene takes the first branch and its
mRNA is never indexed.

Consequences, both measured:

* `get_ref_ids_miniprot` cannot map a miniprot hit on that gene back to a
  reference gene, so every rescue candidate for it is abandoned. 12 genes in
  the rice reference; 0 in human RefSeq.
* `Lifton_feature.children` stays empty, so the childless-gene counter reads
  "the reference had no children either" and does not report a real loss.
  Rice's `gene-OrsajCp001` -- reference children mRNA + 2 exons, emitted as a
  bare gene line -- was missed exactly this way.

Note this is NOT the bacterial `gene -> CDS` shape of GH #37: that has no exon
rows, so it takes the second branch and indexes correctly. Verified directly
against `issue37/bact.gff`.
"""
import tempfile
import types

import pytest

from lifton import annotation, lifton_utils


@pytest.fixture
def organellar_reference(tmp_path):
    """gene -> its own exons, AND gene -> mRNA -> the same exons again."""
    path = tmp_path / "organellar.gff"
    path.write_text(
        "##gff-version 3\n"
        "cp1\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-CP001;gene_biotype=protein_coding\n"
        "cp1\tRefSeq\texon\t100\t199\t.\t+\t.\tID=id-CP001-1;Parent=gene-CP001\n"
        "cp1\tRefSeq\texon\t800\t900\t.\t+\t.\tID=id-CP001-2;Parent=gene-CP001\n"
        "cp1\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-CP001;Parent=gene-CP001\n"
        "cp1\tRefSeq\texon\t100\t199\t.\t+\t.\tID=exon-CP001-1;Parent=rna-CP001\n"
        "cp1\tRefSeq\texon\t800\t900\t.\t+\t.\tID=exon-CP001-2;Parent=rna-CP001\n"
        "cp1\tRefSeq\tCDS\t100\t199\t.\t+\t0\tID=cds-CP001;Parent=rna-CP001\n"
        "cp1\tRefSeq\tCDS\t800\t900\t.\t+\t2\tID=cds-CP001;Parent=rna-CP001\n"
    )
    return path


def _build(path, tmp_path):
    ref = annotation.Annotation(str(path), False, False, "create_unique",
                                None, True, False)
    args = types.SimpleNamespace(evaluation_liftoff_chm13=False, features=None,
                                 gene_only=False, debug=False,
                                 annotation_database="RefSeq")
    out = tempfile.mkdtemp(dir=str(tmp_path))
    returned = lifton_utils.get_ref_liffover_features(["gene"], ref, out, args)
    by_name = {}
    for value in returned:
        if isinstance(value, dict):
            by_name.setdefault(type(value).__name__, []).append(value)
    return returned


class TestOrganellarGeneIsIndexed:
    def test_its_transcript_maps_back_to_its_gene(self, organellar_reference,
                                                  tmp_path):
        """Without this, every miniprot hit on the gene is unrescuable."""
        returned = _build(organellar_reference, tmp_path)
        reverse = next(d for d in returned
                       if isinstance(d, dict) and "rna-CP001" in d
                       ) if any(isinstance(d, dict) and "rna-CP001" in d
                                for d in returned) else None
        assert reverse is not None, (
            "no returned index maps rna-CP001 to its gene; a miniprot hit on "
            "this gene cannot be resolved and is dropped in silence")
        assert reverse["rna-CP001"] == "gene-CP001"

    def test_the_gene_records_that_it_has_children(self, organellar_reference,
                                                   tmp_path):
        """`Lifton_feature.children` is what the childless-gene counter reads."""
        returned = _build(organellar_reference, tmp_path)
        features = next(d for d in returned
                        if isinstance(d, dict) and "gene-CP001" in d
                        and hasattr(next(iter(d.values())), "children"))
        assert features["gene-CP001"].children, (
            "the reference plainly gives this gene an mRNA, so emitting it as "
            "a bare gene line is a loss and must be countable")


class TestThreeLevelIsUnchanged:
    def test_an_ordinary_gene_still_indexes(self, tmp_path):
        path = tmp_path / "ordinary.gff"
        path.write_text(
            "##gff-version 3\n"
            "c1\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-A;gene_biotype=protein_coding\n"
            "c1\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-A;Parent=gene-A\n"
            "c1\tRefSeq\texon\t100\t199\t.\t+\t.\tID=exon-A-1;Parent=rna-A\n"
            "c1\tRefSeq\tCDS\t100\t199\t.\t+\t0\tID=cds-A;Parent=rna-A\n"
        )
        returned = _build(path, tmp_path)
        assert any(isinstance(d, dict) and d.get("rna-A") == "gene-A"
                   for d in returned)
