"""A bare gene line is only worth reporting when the reference had children.

`genes_emitted_without_children` exists because the Liftoff `-copies`
resolution bug emitted ~4,400 genes with no transcript, exon or CDS across the
benchmark corpus and nothing ever added them up. As written it counts every
emitted gene with an empty transcript dict -- which includes the thousands of
single-row pseudogenes RefSeq genuinely declares with no children at all.

Measured on cycle-3 output, comparing each emitted childless gene against the
reference that produced it:

| genome | reported | real |
|---|---:|---:|
| human -> CHM13 | 10,720 | 0 |
| dog -> cat | 272 | 0 |
| human -> zebrafish | 113 | 0 |
| rice | 25 | 2 |

11,130 reported, 2 real. RefSeq alone declares 10,626 childless pseudogenes in
the CHM13 reference. A counter that is 99.98 % false positives cannot show a
recurrence of the bug it was built for, which is the entire point of it.
"""
from types import SimpleNamespace

import pytest

from lifton import lifton_class, locus_pipeline


def _ctx(args, ref_features_dict):
    return SimpleNamespace(args=args, ref_features_dict=ref_features_dict)


def _gene(gene_id, ref_gene_id):
    gene = lifton_class.Lifton_GENE.__new__(lifton_class.Lifton_GENE)
    gene.entry = SimpleNamespace(id=gene_id)
    gene.transcripts = {}
    gene.ref_gene_id = ref_gene_id
    return gene


def _ref(feature_id, child_ids):
    feature = lifton_class.Lifton_feature(feature_id)
    for child in child_ids:
        feature.children.add(child)
    return feature


class TestOnlyRealLossesAreCounted:
    def test_a_reference_childless_pseudogene_is_not_counted(self):
        """RefSeq declares 10,626 of these in the CHM13 reference alone.
        Emitting one with no children is faithful, not a loss."""
        args = SimpleNamespace()
        ref = {"gene-PSEUDO": _ref("gene-PSEUDO", [])}
        locus_pipeline._record_childless_gene(
            _ctx(args, ref), _gene("gene-PSEUDO", "gene-PSEUDO"))
        assert getattr(args, "_childless_gene_count", 0) == 0

    def test_a_gene_whose_reference_had_children_is_counted(self):
        """The `-copies` shape: the reference gave it a transcript and we
        emitted a bare gene line."""
        args = SimpleNamespace()
        ref = {"gene-REAL": _ref("gene-REAL", ["rna-REAL"])}
        locus_pipeline._record_childless_gene(
            _ctx(args, ref), _gene("gene-REAL", "gene-REAL"))
        assert getattr(args, "_childless_gene_count", 0) == 1
        assert getattr(args, "_childless_gene_ids", []) == ["gene-REAL"]

    def test_a_gene_with_children_is_never_counted(self):
        args = SimpleNamespace()
        ref = {"gene-REAL": _ref("gene-REAL", ["rna-REAL"])}
        gene = _gene("gene-REAL", "gene-REAL")
        gene.transcripts = {"rna-REAL": object()}
        locus_pipeline._record_childless_gene(_ctx(args, ref), gene)
        assert getattr(args, "_childless_gene_count", 0) == 0

    def test_an_unknown_reference_gene_is_counted(self):
        """If we cannot tell, say so rather than silently excusing it -- the
        bug this counter exists for is exactly an id that failed to resolve."""
        args = SimpleNamespace()
        locus_pipeline._record_childless_gene(
            _ctx(args, {}), _gene("gene-X_1", "gene-X_1"))
        assert getattr(args, "_childless_gene_count", 0) == 1

    def test_a_copy_resolves_to_its_reference_base(self):
        """`-copies` emits `gene-X_1`; the reference only declares `gene-X`."""
        args = SimpleNamespace()
        ref = {"gene-X": _ref("gene-X", [])}
        locus_pipeline._record_childless_gene(
            _ctx(args, ref), _gene("gene-X_1", "gene-X"))
        assert getattr(args, "_childless_gene_count", 0) == 0


class TestTheRawCountIsStillAvailable:
    def test_every_bare_gene_line_is_tallied_separately(self):
        """The faithful ones are still worth a number, just not that one."""
        args = SimpleNamespace()
        ref = {"gene-PSEUDO": _ref("gene-PSEUDO", []),
               "gene-REAL": _ref("gene-REAL", ["rna-REAL"])}
        ctx = _ctx(args, ref)
        locus_pipeline._record_childless_gene(ctx, _gene("gene-PSEUDO", "gene-PSEUDO"))
        locus_pipeline._record_childless_gene(ctx, _gene("gene-REAL", "gene-REAL"))
        assert getattr(args, "_bare_gene_line_count", 0) == 2
        assert getattr(args, "_childless_gene_count", 0) == 1
