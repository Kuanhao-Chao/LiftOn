"""The genetic code an annotation declares (``transl_table``).

Everything LiftOn does with a coding sequence -- extracting the reference
protein, translating the lifted model, searching for a rescue ORF, completing a
terminal stop -- used the standard code, whatever the annotation said. On a
vertebrate mitochondrial gene (``transl_table=2``) that reads every TGA
tryptophan as a stop, so the reference protein is wrong, the lifted protein is
wrong in the same way, ``stop_codon_gain`` fires, and the ORF search replaces a
correct CDS with the short spurious ORF in front of the first TGA.

Measured on the published CHM13 annotation before this fix: 4 of the 13 human
mitochondrial genes had a truncated CDS (ND2 1042 -> 144 bp, ND4 1378 -> 177,
CYTB 1141 -> 564, COX3 784 -> 405) and the other 9 carried protein identities
from 0.012 to 0.739 where the truth is 1.000.

The table-1 path must not move: an annotation that declares nothing, or
declares 1, has to translate and scan exactly as every prior release did.
"""
from __future__ import annotations

import textwrap
import types

import pytest
from Bio.Seq import Seq

from lifton import coding, lifton, lifton_class, orf_completion
from lifton.exceptions import LiftOnInputError
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


# 13 codons: ATG, 5x GCT, TGA, 5x GCT, TAA. The TGA is a stop under table 1 and
# tryptophan under table 2 -- the whole mitochondrial problem in 39 bases.
MITO_CDS = "ATG" + "GCT" * 5 + "TGA" + "GCT" * 5 + "TAA"
PROTEIN_TABLE_1 = "MAAAAA*AAAAA*"
PROTEIN_TABLE_2 = "MAAAAAWAAAAA*"


def _feature(attributes, featuretype="CDS", identifier="cds1"):
    return types.SimpleNamespace(featuretype=featuretype, id=identifier,
                                 attributes=attributes)


class TestCodonTables:
    def test_stop_codons_follow_the_table(self):
        assert coding.stop_codons(1) == {"TAA", "TAG", "TGA"}
        # TGA is tryptophan here; AGA and AGG are the stops instead.
        assert coding.stop_codons(2) == {"TAA", "TAG", "AGA", "AGG"}
        # Table 11 differs from table 1 only in its start codons, so a plant
        # plastid CDS translates identically -- worth pinning so nobody claims
        # this fix changed those genes.
        assert coding.stop_codons(11) == coding.stop_codons(1)

    def test_default_is_the_standard_code(self):
        assert coding.DEFAULT_TRANSL_TABLE == 1
        assert coding.translate(MITO_CDS) == PROTEIN_TABLE_1

    @pytest.mark.parametrize("table, expected",
                            [(1, PROTEIN_TABLE_1), (2, PROTEIN_TABLE_2),
                             (11, PROTEIN_TABLE_1)])
    def test_translate_matches_biopython(self, table, expected):
        assert coding.translate(MITO_CDS, table) == expected
        assert coding.translate(MITO_CDS, table) == str(
            Seq(MITO_CDS).translate(table=table))

    def test_translate_drops_a_trailing_partial_codon(self):
        # The established LiftOn result, made explicit so it cannot drift.
        assert coding.translate(MITO_CDS + "AT", 2) == PROTEIN_TABLE_2
        assert coding.translate("AT", 1) == ""
        assert coding.translate("", 1) == ""

    def test_reverse_strand_is_the_callers_business(self):
        # Sequences reach translate() already in transcript orientation.
        revcomp = str(Seq(MITO_CDS).reverse_complement())
        assert coding.translate(str(Seq(revcomp).reverse_complement()), 2) == \
            PROTEIN_TABLE_2

    def test_an_unknown_table_is_rejected_by_name(self):
        with pytest.raises(LiftOnInputError) as error:
            coding.parse_transl_table("99", "cds-XYZ")
        assert "99" in str(error.value) and "cds-XYZ" in str(error.value)
        with pytest.raises(LiftOnInputError):
            coding.parse_transl_table("mitochondrial")

    def test_values_arrive_as_gffutils_lists(self):
        assert coding.table_from_attributes({"transl_table": ["2"]}) == 2
        assert coding.table_from_attributes({"transl_table": "2"}) == 2
        assert coding.table_from_attributes({"transl_table": []}) is None
        assert coding.table_from_attributes({"gene": ["ND1"]}) is None
        assert coding.table_from_attributes(None) is None


class TestResolveTranslTable:
    def test_declared_on_the_cds(self):
        features = [_feature({"transl_table": ["2"]})]
        assert coding.resolve_transl_table(features) == 2

    def test_nothing_declared_means_the_standard_code(self):
        assert coding.resolve_transl_table([_feature({"gene": ["x"]})]) == 1
        assert coding.resolve_transl_table([]) == 1
        assert coding.resolve_transl_table(None) == 1

    def test_only_cds_rows_are_consulted(self):
        # An exon carrying the attribute is not the transcript's coding code.
        exon = _feature({"transl_table": ["2"]}, featuretype="exon")
        assert coding.resolve_transl_table([exon]) == 1

    def test_first_in_transcript_order_wins_and_the_conflict_is_reported(self, capsys):
        features = [_feature({"transl_table": ["2"]}, identifier="a"),
                    _feature({"transl_table": ["11"]}, identifier="b")]
        assert coding.resolve_transl_table(features, context="rna-X") == 2
        assert "more than one" in capsys.readouterr().err


class TestTranscriptCarriesItsTable:
    def _transcript(self, cds_attributes=None):
        entry = types.SimpleNamespace(
            attributes={"ID": ["rna-X"]}, id="rna-X", source="test",
            seqid="chrM", start=1, end=39, strand="+")
        trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
        trans.entry, trans.exons, trans.exon_dic = entry, [], {}
        if cds_attributes is not None:
            trans._transl_table = coding.table_from_attributes(cds_attributes)
        return trans

    def test_a_model_with_no_declared_code_uses_the_standard_one(self):
        assert self._transcript().transl_table() == 1
        assert self._transcript({}).transl_table() == 1

    def test_a_declared_code_reaches_the_translation(self):
        trans = self._transcript({"transl_table": ["2"]})
        assert trans.transl_table() == 2
        assert trans.translate_coding_seq(MITO_CDS) == PROTEIN_TABLE_2
        assert self._transcript().translate_coding_seq(MITO_CDS) == PROTEIN_TABLE_1

    def test_add_cds_captures_the_code_from_the_first_reference_cds(self):
        trans = self._transcript()
        cds = types.SimpleNamespace(
            featuretype="CDS", id="cds1", seqid="chrM", start=1, end=39,
            strand="+", frame="0", attributes={"transl_table": ["2"],
                                               "Parent": ["rna-X"]})
        trans.add_exon(types.SimpleNamespace(
            featuretype="exon", id="exon1", seqid="chrM", start=1, end=39,
            strand="+", frame=".", attributes={"Parent": ["rna-X"]}))
        trans.add_cds(cds)
        assert trans.transl_table() == 2


class TestStopCompletionFollowsTheTable:
    def test_the_stop_set_comes_from_the_model(self):
        standard = types.SimpleNamespace(transl_table=lambda: 1)
        mito = types.SimpleNamespace(transl_table=lambda: 2)
        assert "TGA" in orf_completion._stop_codons(standard)
        assert "TGA" not in orf_completion._stop_codons(mito)
        assert "AGA" in orf_completion._stop_codons(mito)

    def test_a_model_without_the_method_keeps_the_standard_stops(self):
        assert orf_completion._stop_codons(object()) == orf_completion.STOP_CODONS


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


def _build_workspace(work, *, declare_table):
    """A single-exon mitochondrial-shaped gene, lifted onto an identical
    target. The DNA lift is exact, so anything but a perfect model is the
    genetic code being read wrong."""
    work.mkdir(parents=True, exist_ok=True)
    chromosome = ["A"] * 600
    chromosome[100:100 + len(MITO_CDS)] = MITO_CDS
    sequence = ">chrM\n" + _wrap("".join(chromosome))
    (work / "ref.fa").write_text(sequence)
    (work / "tgt.fa").write_text(sequence)
    table = ";transl_table=2" if declare_table else ""
    rows = (
        "chrM\t{src}\tgene\t101\t139\t.\t+\t.\tID=gene-ND;gene_biotype=protein_coding\n"
        "chrM\t{src}\tmRNA\t101\t139\t.\t+\t.\tID=rna-ND;Parent=gene-ND\n"
        "chrM\t{src}\texon\t101\t139\t.\t+\t.\tID=exon-ND;Parent=rna-ND\n"
        "chrM\t{src}\tCDS\t101\t139\t.\t+\t0\tID=cds-ND;Parent=rna-ND" + table + "\n")
    (work / "ref.gff3").write_text("##gff-version 3\n" + rows.format(src="test"))
    (work / "liftoff.gff3").write_text(
        "##gff-version 3\n" + rows.format(src="Liftoff"))
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chrM\tminiprot\tmRNA\t101\t139\t.\t+\t.\tID=MP1;Rank=1;Identity=1.0000;"
        "Target=rna-ND 1 13\n"
        "chrM\tminiprot\tCDS\t101\t139\t.\t+\t0\tID=MP1.c1;Parent=MP1\n")
    (work / "out").mkdir()
    return work


def _run(work):
    output = work / "out" / "lifton.gff3"
    lifton.run_all_lifton_steps(lifton.parse_args([
        str(work / "tgt.fa"), str(work / "ref.fa"),
        "-g", str(work / "ref.gff3"), "-L", str(work / "liftoff.gff3"),
        "-M", str(work / "miniprot.gff3"), "-o", str(output),
        "-ad", "RefSeq", "--force"]))
    return work, output.read_text()


def _mrna_attributes(text):
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "mRNA":
            return dict(kv.split("=", 1) for kv in fields[8].split(";")
                        if "=" in kv)
    raise AssertionError("no mRNA row in the output")


def _cds_span(text):
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "CDS":
            return int(fields[3]), int(fields[4])
    raise AssertionError("no CDS row in the output")


class TestMitochondrialGeneEndToEnd:
    def test_the_declared_table_keeps_the_cds_and_scores_it_correctly(
            self, tmp_path, hermetic_pipeline):
        work, text = _run(_build_workspace(tmp_path / "mito", declare_table=True))
        attributes = _mrna_attributes(text)
        assert attributes["protein_identity"] == "1.000"
        assert "stop_codon_gain" not in attributes.get("mutation", "")
        assert _cds_span(text) == (101, 139)
        proteins = (work / "out" / "lifton_output" / "intermediate_files"
                    / "proteins.fa").read_text()
        assert PROTEIN_TABLE_2.rstrip("*") in proteins.replace("\n", "")

    def test_without_the_attribute_the_standard_code_still_applies(
            self, tmp_path, hermetic_pipeline):
        # The pre-fix behaviour, kept exactly: this is what a table-1
        # annotation (and every annotation that declares nothing) must do.
        work, text = _run(_build_workspace(tmp_path / "std", declare_table=False))
        attributes = _mrna_attributes(text)
        assert float(attributes["protein_identity"]) < 1.0
        proteins = (work / "out" / "lifton_output" / "intermediate_files"
                    / "proteins.fa").read_text()
        assert PROTEIN_TABLE_1.rstrip("*") in proteins.replace("\n", "")


class TestMiniprotSpeaksOneCodePerInvocation:
    """miniprot takes one ``-T`` per run, so a reference that declares more
    than one genetic code needs one invocation per code.

    Measured on the 13 human mitochondrial proteins against CHM13 chrM: run
    under the standard code instead of the table 2 they declare, every one of
    them loses identity (mean 1.0000 -> 0.9302; COX1 picks up 16 spurious
    in-frame stops) and 4 of the 13 hits land on different coordinates.
    """

    def _proteins(self, tmp_path, tables):
        from lifton import extract_sequence
        path = tmp_path / "proteins.fa"
        path.write_text("".join(f">{name}\nMAAA\n" for name in tables))
        extract_sequence.write_transl_table_sidecar(
            str(path), {name: table for name, table in tables.items()
                        if table != 1})
        return str(path)

    def test_a_single_code_reference_is_untouched(self, tmp_path):
        from lifton import run_miniprot
        path = self._proteins(tmp_path, {"a": 1, "b": 1})
        assert run_miniprot.transl_table_groups(path) == [(1, path)]
        # No sidecar at all -- a user-supplied -P file -- behaves the same.
        plain = tmp_path / "user.faa"
        plain.write_text(">a\nMAAA\n")
        assert run_miniprot.transl_table_groups(str(plain)) == [(1, str(plain))]

    def test_groups_are_split_with_the_standard_code_first(self, tmp_path):
        from lifton import run_miniprot
        path = self._proteins(tmp_path, {"nuc1": 1, "mito": 2, "nuc2": 1,
                                         "plastid": 11})
        groups = run_miniprot.transl_table_groups(path)
        assert [table for table, _ in groups] == [1, 2, 11]
        contents = {table: open(fasta).read() for table, fasta in groups}
        assert ">nuc1" in contents[1] and ">nuc2" in contents[1]
        assert contents[2].startswith(">mito") and ">nuc1" not in contents[2]
        assert contents[11].startswith(">plastid")
        # Every protein is placed exactly once.
        assert sum(text.count(">") for text in contents.values()) == 4

    def test_the_standard_group_command_is_unchanged(self):
        from lifton import run_miniprot
        base = run_miniprot._build_miniprot_command(
            "miniprot", "tgt.fa", "p.faa", "", 1)
        assert base == run_miniprot._build_miniprot_command(
            "miniprot", "tgt.fa", "p.faa", "", 1, table=1)
        assert "-T" not in base and "-P" not in base

    def test_a_non_standard_group_gets_its_table_and_its_own_ids(self):
        from lifton import run_miniprot
        command = run_miniprot._build_miniprot_command(
            "miniprot", "tgt.fa", "p.table2.faa", "", 1, table=2)
        assert command[command.index("-T") + 1] == "2"
        # Both invocations number their hits from MP000001, so the second
        # group needs its own prefix or the ids collide.
        assert command[command.index("-P") + 1] == "MPT2"

    def test_explicit_user_options_win(self):
        from lifton import run_miniprot
        command = run_miniprot._build_miniprot_command(
            "miniprot", "tgt.fa", "p.faa", "-T 5 -P X", 1, table=2)
        assert command.count("-T") == 1 and command.count("-P") == 1
        assert "MPT2" not in command

    def test_the_sidecar_is_removed_when_nothing_is_non_standard(self, tmp_path):
        from lifton import extract_sequence
        path = str(tmp_path / "proteins.fa")
        open(path, "w").write(">a\nM\n")
        extract_sequence.write_transl_table_sidecar(path, {"a": 2})
        assert extract_sequence.read_transl_table_sidecar(path) == {"a": 2}
        # A stale sidecar from an earlier run must not steer a later one.
        extract_sequence.write_transl_table_sidecar(path, {})
        assert extract_sequence.read_transl_table_sidecar(path) == {}


class TestMiniprotDerivedModels:
    """A miniprot row carries no transl_table, so every model built from one --
    Step 8, the candidate scaffold, the rescue and its isoforms, cross-locus --
    was translated, ORF-scanned and stop-completed with the standard code even
    when its reference transcript declares 2 or 5: TGA read as a stop."""

    def teardown_method(self):
        from lifton import coding
        coding.install_transcript_tables({})

    def _miniprot_model(self, ref_id):
        from lifton import lifton_class
        trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
        trans.ref_tran_id = ref_id
        return trans

    def test_the_reference_transcript_code_is_used(self):
        from lifton import coding, orf_completion
        coding.install_transcript_tables({"rna-ND4": 2})
        trans = self._miniprot_model("rna-ND4")
        assert trans.transl_table() == 2
        assert trans.translate_coding_seq("ATGTGAAGA") == "MW*"
        assert "AGA" in orf_completion._stop_codons(trans)

    def test_a_transcript_that_declares_nothing_keeps_the_standard_code(self):
        from lifton import coding
        coding.install_transcript_tables({"rna-ND4": 2})
        trans = self._miniprot_model("rna-OTHER")
        assert trans.transl_table() == coding.DEFAULT_TRANSL_TABLE
        assert trans.translate_coding_seq("ATGTGAAGA") == "M*R"
