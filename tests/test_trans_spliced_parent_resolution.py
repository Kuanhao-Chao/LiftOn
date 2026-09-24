"""A repeated trans-spliced gene ID must not bind children across seqids."""

import io
from collections import defaultdict
from types import SimpleNamespace

import gffutils

from lifton import coreutils, lifton_class
from lifton.gff3_validator import validate_gff3_file
from lifton.locus_pipeline import HierarchyBatchLoader, _ThreadLocalCtxFactory
from lifton.parent_resolution import bind_same_seqid_parents


def _database(tmp_path, *, ambiguous=False):
    lines = [
        "##gff-version 3",
        "chrA\tLiftoff\tgene\t10\t30\t.\t+\t.\tID=g;part=3;exception=trans-splicing",
        "chrB\tLiftoff\tgene\t100\t200\t.\t+\t.\tID=g;part=1;exception=trans-splicing",
    ]
    if ambiguous:
        lines.append("chrB\tLiftoff\tgene\t100\t200\t.\t+\t.\tID=g;part=2;exception=trans-splicing")
    lines += [
        "chrB\tLiftoff\tmRNA\t120\t180\t.\t+\t.\tID=t;Parent=g;exception=trans-splicing",
        "chrB\tLiftoff\texon\t120\t180\t.\t+\t.\tID=e;Parent=t",
        "chrB\tLiftoff\tCDS\t120\t180\t.\t+\t0\tID=c;Parent=t",
    ]
    path = tmp_path / "lift.gff3"
    path.write_text("\n".join(lines) + "\n")
    return gffutils.create_db(
        str(path), dbfn=str(tmp_path / "lift.db"), force=True,
        merge_strategy="create_unique", disable_infer_genes=True,
        disable_infer_transcripts=True,
    )


def test_duplicate_id_transcript_moves_to_unique_same_seqid_parent(tmp_path):
    database = _database(tmp_path)
    roots = list(database.features_of_type("gene"))
    first, second = roots
    assert first.seqid == "chrA" and second.seqid == "chrB"
    assert [child.id for child in database.children(first, level=1)] == ["t"]
    assert list(database.children(second, level=1)) == []

    repaired = bind_same_seqid_parents(database, ["gene"])
    assert repaired.repaired_branches == 1
    assert list(repaired.children(first, level=1)) == []
    assert [child.id for child in repaired.children(second, level=1)] == ["t"]
    assert [child.id for child in repaired.children(second, featuretype="CDS")] == ["c"]
    assert [child.id for child in repaired.children("g_1", featuretype="exon")] == ["e"]
    assert {child.id for child in repaired.children("t", level=1)} == {"e", "c"}


def test_batched_locus_loader_sees_corrected_children(tmp_path):
    database = _database(tmp_path)

    class BatchedDatabase:
        def __getitem__(self, key):
            return database[key]

        def features_of_type(self, *args, **kwargs):
            return database.features_of_type(*args, **kwargs)

        def children(self, *args, **kwargs):
            return database.children(*args, **kwargs)

        def children_batched_features(self, anchors, **kwargs):
            return {anchor: tuple(database.children(anchor, **kwargs)) for anchor in anchors}

    repaired = bind_same_seqid_parents(BatchedDatabase(), ["gene"])
    loader = HierarchyBatchLoader(repaired, slots=2)
    assert loader.batched
    children = loader.children_many(["g", "g_1"], level=1)
    assert children["g"] == ()
    assert [child.id for child in children["g_1"]] == ["t"]


def test_thread_reopen_preserves_parent_repair(tmp_path):
    database = _database(tmp_path)
    repaired = bind_same_seqid_parents(database, ["gene"])
    worker = _ThreadLocalCtxFactory._open_thread_db(repaired, database.dbfn)
    try:
        assert worker is not repaired
        assert worker.conn is not database.conn
        assert worker.repaired_branches == 1
        assert list(worker.children("g", level=1)) == []
        assert [child.id for child in worker.children("g_1", level=1)] == ["t"]
    finally:
        _ThreadLocalCtxFactory._close_db(worker)


def test_trans_spliced_part_attribute_follows_lifted_fragment(tmp_path):
    database = _database(tmp_path)
    first = database["g"]
    second = database["g_1"]

    class Journal:
        def allocate_gene_copy(self, *args):
            return 1

    gene = lifton_class.Lifton_GENE(
        "g", second, dict(first.attributes), {}, {},
        SimpleNamespace(annotation_database="RefSeq"), state_journal=Journal(),
    )
    assert gene.entry.attributes["part"] == ["1"]
    assert gene.entry.id == "g_1"


def test_rebound_family_serializes_without_cross_seqid_parent(tmp_path):
    database = _database(tmp_path)
    repaired = bind_same_seqid_parents(database, ["gene"])
    first, second = list(database.features_of_type("gene"))

    class Journal:
        def __init__(self, copy_num):
            self.copy_num = copy_num

        def allocate_gene_copy(self, *args):
            return self.copy_num

    output = io.StringIO()
    stats = defaultdict(dict)
    for root, copy_num in ((first, 0), (second, 1)):
        gene = lifton_class.Lifton_GENE(
            "g", coreutils.clone_feature(root), dict(root.attributes), {}, {},
            SimpleNamespace(annotation_database="RefSeq"),
            state_journal=Journal(copy_num),
        )
        for child in repaired.children(root, level=1):
            trans = gene.add_transcript("t", coreutils.clone_feature(child), dict(child.attributes))
            for exon in repaired.children(child, featuretype="exon", level=1):
                gene.add_exon(trans.entry.id, coreutils.clone_feature(exon))
            for cds in repaired.children(child, featuretype="CDS", level=1):
                gene.add_cds(trans.entry.id, coreutils.clone_feature(cds))
        assert gene.write_entry(output, stats)

    path = tmp_path / "rebound.gff3"
    path.write_text("##gff-version 3\n" + output.getvalue())
    result = validate_gff3_file(path)
    assert result.errors == []
    rows = [line for line in output.getvalue().splitlines() if "\tmRNA\t" in line]
    assert len(rows) == 1 and rows[0].startswith("chrB\t")
    assert "Parent=g_1" in rows[0]


def test_ambiguous_same_seqid_family_is_left_as_bound_and_counted(tmp_path):
    """Two fragments could hold the transcript, so neither is chosen. The run
    must go on: raising here aborted a whole genome over one gene that every
    earlier release lifted (with its cross-sequence parent)."""
    database = _database(tmp_path, ambiguous=True)
    counts = {}
    bound = bind_same_seqid_parents(database, ["gene"], counts)
    assert bound is database
    assert counts == {"repaired": 0, "ambiguous": 1}
    assert [child.id for child in bound.children("g", level=1)] == ["t"]


def test_counts_report_a_repair(tmp_path):
    counts = {}
    bound = bind_same_seqid_parents(_database(tmp_path), ["gene"], counts)
    assert counts == {"repaired": 1, "ambiguous": 0}
    assert bound.repaired_branches == 1


def test_unique_ids_keep_original_database_object(tmp_path):
    database = _database(tmp_path)
    assert bind_same_seqid_parents(database, ["mRNA"]) is database


def _family_on_one_sequence(tmp_path, first_root_span):
    """drosophila mod(mdg4) in miniature: fragments of one trans-spliced gene
    on the SAME sequence, the transcript lying inside the second."""
    start, end = first_root_span
    lines = [
        "##gff-version 3",
        f"chrA\tLiftoff\tgene\t{start}\t{end}\t.\t-\t.\tID=g;part=2;exception=trans-splicing",
        "chrA\tLiftoff\tgene\t100\t200\t.\t-\t.\tID=g;part=1;exception=trans-splicing",
        "chrA\tLiftoff\tmRNA\t120\t180\t.\t-\t.\tID=t;Parent=g;exception=trans-splicing",
        "chrA\tLiftoff\texon\t120\t180\t.\t-\t.\tID=e;Parent=t",
        "chrA\tLiftoff\tCDS\t120\t180\t.\t-\t0\tID=c;Parent=t",
    ]
    path = tmp_path / "same_sequence.gff3"
    path.write_text("\n".join(lines) + "\n")
    return gffutils.create_db(
        str(path), dbfn=str(tmp_path / "same_sequence.db"), force=True,
        merge_strategy="create_unique", disable_infer_genes=True,
        disable_infer_transcripts=True,
    )


def test_same_sequence_transcript_moves_to_the_fragment_that_contains_it(tmp_path):
    """Bound to a fragment that does not contain it, the transcript made that
    fragment's span normalise to the other fragment's coordinates -- written
    twice at one locus, and no gene row at its own."""
    database = _family_on_one_sequence(tmp_path, (10, 30))
    counts = {}
    bound = bind_same_seqid_parents(database, ["gene"], counts)
    assert counts == {"repaired": 1, "ambiguous": 0}
    assert list(bound.children("g", level=1)) == []
    assert [child.id for child in bound.children("g_1", level=1)] == ["t"]
    assert {child.id for child in bound.children("g_1")} == {"t", "e", "c"}


def test_a_transcript_its_bound_fragment_contains_never_moves(tmp_path):
    """Overlapping fragments can both contain a child. Its bound root holding
    it is the database's own answer, and nothing here second-guesses it."""
    database = _family_on_one_sequence(tmp_path, (90, 210))
    counts = {}
    bound = bind_same_seqid_parents(database, ["gene"], counts)
    assert counts == {"repaired": 0, "ambiguous": 0}
    assert bound is database


def test_a_declared_id_missing_from_the_database_is_skipped():
    """gffutils raises FeatureNotFoundError -- not KeyError -- for a missing
    key, so the overlay's guard let it escape and abort the lift."""
    from types import SimpleNamespace
    import gffutils
    from lifton import parent_resolution

    class Database:
        def features_of_type(self, featuretype):
            return iter([SimpleNamespace(id="gene-A_1", attributes={"ID": ["gene-A"]},
                                         featuretype="gene", seqid="chr1", start=1,
                                         end=10, strand="+")])

        def __getitem__(self, key):
            raise gffutils.exceptions.FeatureNotFoundError(key)

    overlay = parent_resolution.SameSeqidParentOverlay(Database(), ["gene"])
    assert overlay.ambiguous == []
