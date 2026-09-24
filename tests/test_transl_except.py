"""transl_except: codons the reference declares to translate differently.

Reported by a user lifting MANE to CHM13: SEPHS2 (a selenoprotein; UGA recoded
as selenocysteine, declared with ``transl_except=(pos:complement(...),aa:Sec)``)
came out with its N-terminus cut off. Measured on the GRCh38 -> CHM13 RefSeq
lift: all 25 human selenoprotein genes were affected, 53 transcripts at a mean
protein identity of 0.666 on a same-species lift, and every carried
transl_except kept the reference's coordinates.

The reference protein and the lifted protein both carry ``*`` at the Sec codon.
The identity cut-off read it as a premature stop (an identical selenoprotein
scored residue/length), find_variants called it stop_codon_gain, and the
correct Liftoff model lost to a miniprot fragment. These tests pin the
read-through, the placement arithmetic, and the output remap.
"""

from __future__ import annotations

import json
import textwrap
from types import SimpleNamespace

import pytest

from lifton import coding, get_id_fraction, transl_except


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

class TestParse:
    def test_decoded_values_one_per_exception(self):
        # What gffutils returns for SELENON (two Sec).
        values = ["(pos:25802093..25802095,aa:Sec)", "(pos:25812789..25812791,aa:Sec)"]
        assert coding.parse_transl_except(values) == [
            ("25802093..25802095", "Sec"), ("25812789..25812791", "Sec")]

    def test_raw_encoded_value(self):
        raw = ["(pos:complement(30445548..30445550)%2Caa:Sec)"]
        assert coding.parse_transl_except(raw) == [
            ("complement(30445548..30445550)", "Sec")]

    def test_split_fragments_are_rejoined(self):
        # A parser that splits on literal commas cuts one exception in two.
        assert coding.parse_transl_except(
            ["(pos:complement(join(559..560", "661..661))", "aa:Sec)"]) == [
            ("complement(join(559..560,661..661))", "Sec")]

    @pytest.mark.parametrize("location,positions", [
        ("140..142", [140, 141, 142]),
        ("complement(140..142)", [142, 141, 140]),
        ("join(10..11,20..20)", [10, 11, 20]),
        ("complement(join(559..560,661..661))", [661, 560, 559]),
        ("join(complement(661..661),complement(559..560))", [661, 560, 559]),
        ("10404", [10404]),
        ("<5..>7", [5, 6, 7]),
    ])
    def test_locations(self, location, positions):
        assert coding.parse_location(location) == positions

    @pytest.mark.parametrize("bad", ["abc", "9..3", "join(1..2"])
    def test_malformed_locations_raise(self, bad):
        with pytest.raises(ValueError):
            coding.parse_location(bad)

    def test_a_malformed_value_is_skipped_not_fatal(self):
        assert coding.transl_excepts_for(
            ["(pos:bad,aa:Sec)"], [(1, 9)], "+", 0, "MA*", "tx") == ()


# ---------------------------------------------------------------------------
# Placement on the reference protein
# ---------------------------------------------------------------------------

class TestPlacement:
    def test_plus_strand(self):
        # CDS 101-160 + 261-278, Sec codon 107-109: residue 2.
        placed = coding.transl_excepts_for(
            ["(pos:107..109,aa:Sec)"], [(101, 160), (261, 278)], "+", 0,
            "MA*" + "A" * 22 + "*")
        assert placed == (coding.TranslExcept(2, "Sec", 3, True),)

    def test_minus_strand_codon_split_across_an_exon_junction(self):
        placed = coding.transl_excepts_for(
            ["(pos:complement(join(559..560,661..661)),aa:Sec)"],
            [(520, 560), (661, 700)], "-", 0,
            "M" + "A" * 12 + "*" + "A" * 12 + "*")
        assert placed == (coding.TranslExcept(13, "Sec", 3, True),)

    def test_initial_phase_shifts_the_frame(self):
        # With phase 1 the first base is skipped: offset 7 is residue 2.
        placed = coding.transl_excepts_for(
            ["(pos:108..110,aa:Sec)"], [(101, 200)], "+", 1, "MA*" + "A" * 20)
        assert placed[0].residue == 2

    def test_out_of_frame_or_outside_the_cds_is_unplaced(self):
        out_of_frame = coding.transl_excepts_for(
            ["(pos:108..110,aa:Sec)"], [(101, 200)], "+", 0, "M" * 33)
        outside = coding.transl_excepts_for(
            ["(pos:50..52,aa:Sec)"], [(101, 200)], "+", 0, "M" * 33)
        assert out_of_frame[0].residue is None and not out_of_frame[0].readthrough
        assert outside[0].residue is None

    def test_readthrough_only_over_an_internal_stop(self):
        # Met over a non-AUG start, a partial TERM, a Sec at the terminal codon:
        # none of them is read through.
        protein = "LAA*"
        met, term, terminal = coding.transl_excepts_for(
            ["(pos:1..3,aa:Met)", "(pos:13,aa:TERM)", "(pos:10..12,aa:Sec)"],
            [(1, 13)], "+", 0, protein)
        assert (met.residue, met.readthrough) == (0, False)
        assert (term.length, term.readthrough) == (1, False)
        assert (terminal.residue, terminal.readthrough) == (3, False)


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

class TestScoring:
    def test_no_columns_is_the_legacy_score(self):
        for reference, target in [("MA*AA*", "MA*AA*"), ("MAKLV*", "MAK-V*"),
                                  ("MA-AA*", "MAKAA*")]:
            assert (get_id_fraction.get_AA_id_fraction(reference, target)
                    == get_id_fraction.get_AA_id_fraction(reference, target, frozenset()))

    def test_an_identical_selenoprotein_scores_one(self):
        columns = get_id_fraction.readthrough_columns("MA*AA*", "MA*AA*", frozenset({2}))
        assert get_id_fraction.get_AA_id_fraction("MA*AA*", "MA*AA*") == (3, 6)
        assert get_id_fraction.get_AA_id_fraction("MA*AA*", "MA*AA*", columns) == (6, 6)

    def test_an_ncbi_u_in_the_reference_matches_the_lifted_stop(self):
        columns = get_id_fraction.readthrough_columns("MAUAA", "MA*AA", frozenset({2}))
        assert get_id_fraction.get_AA_id_fraction("MAUAA", "MA*AA", columns) == (5, 5)

    def test_a_stop_after_the_declared_one_still_ends_the_scan(self):
        reference, target = "MA*AA*AA", "MA*A*AAA"
        columns = get_id_fraction.readthrough_columns(reference, target, frozenset({2}))
        assert get_id_fraction.get_AA_id_fraction(reference, target, columns) == (4, 8)

    def test_a_model_that_ends_at_the_declared_codon_stops_there(self):
        assert get_id_fraction.readthrough_columns(
            "MA*AA*", "MA*---", frozenset({2})) == frozenset()

    def test_partial_windows(self):
        reference, target = "MA*AA*", "MA*AA*"
        columns = get_id_fraction.readthrough_columns(reference, target, frozenset({2}))
        assert get_id_fraction.get_partial_id_fraction(reference, target, 0, 5) == (3, 5)
        assert get_id_fraction.get_partial_id_fraction(
            reference, target, 0, 5, columns) == (5, 5)

    def test_the_masked_protein_splits_only_on_real_stops(self):
        columns = get_id_fraction.readthrough_columns("MA*AA*", "MA*AA*", frozenset({2}))
        masked = get_id_fraction.mask_readthrough("MA*AA*", "MA*AA*", columns)
        assert masked.split("*") == ["MAUAA", ""]


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

class TestRegistry:
    def teardown_method(self):
        transl_except.clear()

    def test_unknown_transcripts_get_none(self):
        transl_except.install({})
        assert transl_except.readthrough("tx") is None
        assert transl_except.get("tx") is None

    def test_a_residue_the_loaded_protein_disagrees_with_is_not_read_through(self):
        proteins = {"tx": "MAKAA*"}
        transl_except.install(
            {"tx": (coding.TranslExcept(2, "Sec", 3, True),)}, _Proteins(proteins))
        assert transl_except.readthrough("tx") is None

    def test_an_ncbi_u_is_accepted(self):
        transl_except.install(
            {"tx": (coding.TranslExcept(2, "Sec", 3, True),)}, _Proteins({"tx": "MAUAA"}))
        assert transl_except.readthrough("tx") == frozenset({2})


class _Proteins(dict):
    """The slice of the pyfaidx interface the registry uses."""


# ---------------------------------------------------------------------------
# End to end, through the real pipeline
# ---------------------------------------------------------------------------

def _wrap(seq):
    return "\n".join(textwrap.wrap(seq, 60)) + "\n"


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


#: The SEPHS2 shape: selenocysteine early (residue 2 of 26), so a score that
#: stops at it rates the correct model 3/26 and a miniprot fragment wins.
EARLY = "ATG" + "GCT" + "TGA" + "GCT" * 22 + "TAA"
#: Selenocysteine at residue 13 of 27.
MID = "ATG" + "GCT" * 12 + "TGA" + "GCT" * 12 + "TAA"
PLAIN = "ATG" + "GCT" * 25 + "TAA"


def _intron(length):
    return "GT" + "C" * (length - 4) + "AG"


def _genomes():
    ref = ["C"] * 1500
    tgt = ["C"] * 3000

    def put(chrom, start, seq):
        chrom[start - 1:start - 1 + len(seq)] = list(seq)

    # A (+): exon 101-160, intron 100 bp, exon 261-278 | target +1000, intron 150
    put(ref, 101, EARLY[:60] + _intron(100) + EARLY[60:])
    put(tgt, 1101, EARLY[:60] + _intron(150) + EARLY[60:])
    # B (-): Sec codon split across the junction. Transcript order: exon 1
    # (40 nt) + intron + exon 2 (41 nt); genome holds the reverse complement.
    put(ref, 520, _revcomp(MID[:40] + _intron(100) + MID[40:]))
    put(tgt, 1520, _revcomp(MID[:40] + _intron(150) + MID[40:]))
    # C (+): the target has TGC (Cys) where the reference has TGA (Sec).
    put(ref, 801, MID)
    put(tgt, 2001, MID[:39] + "TGC" + MID[42:])
    # D (+): no exception at all.
    put(ref, 1001, PLAIN)
    put(tgt, 2201, PLAIN)
    return "".join(ref), "".join(tgt)


def _gff(source, shift_a, shift_b_exon2, shift_b_exon1, c, d, with_exceptions):
    """One annotation of the four genes. transl_except always carries the
    REFERENCE coordinates, exactly as Liftoff copies it."""
    te_a = ";transl_except=(pos:107..109%2Caa:Sec)" if with_exceptions else ""
    te_b = (";transl_except=(pos:complement(join(559..560%2C661..661))%2Caa:Sec)"
            if with_exceptions else "")
    te_c = ";transl_except=(pos:840..842%2Caa:Sec)" if with_exceptions else ""
    a1, a2 = (101 + shift_a, 160 + shift_a), (261 + shift_a + (50 if shift_a else 0),
                                              278 + shift_a + (50 if shift_a else 0))
    b2 = (520 + shift_b_exon2, 560 + shift_b_exon2)
    b1 = (661 + shift_b_exon1, 700 + shift_b_exon1)
    rows = [
        f"chr1\t{source}\tgene\t{a1[0]}\t{a2[1]}\t.\t+\t.\tID=geneA;gene_biotype=protein_coding",
        f"chr1\t{source}\tmRNA\t{a1[0]}\t{a2[1]}\t.\t+\t.\tID=txA;Parent=geneA",
        f"chr1\t{source}\texon\t{a1[0]}\t{a1[1]}\t.\t+\t.\tID=exonA1;Parent=txA",
        f"chr1\t{source}\texon\t{a2[0]}\t{a2[1]}\t.\t+\t.\tID=exonA2;Parent=txA",
        f"chr1\t{source}\tCDS\t{a1[0]}\t{a1[1]}\t.\t+\t0\tID=cdsA;Parent=txA{te_a}",
        f"chr1\t{source}\tCDS\t{a2[0]}\t{a2[1]}\t.\t+\t0\tID=cdsA;Parent=txA{te_a}",
        f"chr1\t{source}\tgene\t{b2[0]}\t{b1[1]}\t.\t-\t.\tID=geneB;gene_biotype=protein_coding",
        f"chr1\t{source}\tmRNA\t{b2[0]}\t{b1[1]}\t.\t-\t.\tID=txB;Parent=geneB",
        f"chr1\t{source}\texon\t{b2[0]}\t{b2[1]}\t.\t-\t.\tID=exonB2;Parent=txB",
        f"chr1\t{source}\texon\t{b1[0]}\t{b1[1]}\t.\t-\t.\tID=exonB1;Parent=txB",
        f"chr1\t{source}\tCDS\t{b2[0]}\t{b2[1]}\t.\t-\t2\tID=cdsB;Parent=txB{te_b}",
        f"chr1\t{source}\tCDS\t{b1[0]}\t{b1[1]}\t.\t-\t0\tID=cdsB;Parent=txB{te_b}",
        f"chr1\t{source}\tgene\t{c}\t{c + 80}\t.\t+\t.\tID=geneC;gene_biotype=protein_coding",
        f"chr1\t{source}\tmRNA\t{c}\t{c + 80}\t.\t+\t.\tID=txC;Parent=geneC",
        f"chr1\t{source}\texon\t{c}\t{c + 80}\t.\t+\t.\tID=exonC;Parent=txC",
        f"chr1\t{source}\tCDS\t{c}\t{c + 80}\t.\t+\t0\tID=cdsC;Parent=txC{te_c}",
        f"chr1\t{source}\tgene\t{d}\t{d + 80}\t.\t+\t.\tID=geneD;gene_biotype=protein_coding",
        f"chr1\t{source}\tmRNA\t{d}\t{d + 80}\t.\t+\t.\tID=txD;Parent=geneD",
        f"chr1\t{source}\texon\t{d}\t{d + 80}\t.\t+\t.\tID=exonD;Parent=txD",
        f"chr1\t{source}\tCDS\t{d}\t{d + 80}\t.\t+\t0\tID=cdsD;Parent=txD",
    ]
    return "##gff-version 3\n" + "\n".join(rows) + "\n"


def _workspace(root, with_exceptions=True):
    root.mkdir(parents=True)
    ref, tgt = _genomes()
    (root / "ref.fa").write_text(">chr1\n" + _wrap(ref))
    (root / "tgt.fa").write_text(">chr1\n" + _wrap(tgt))
    (root / "ref.gff3").write_text(
        _gff("test", 0, 0, 0, 801, 1001, with_exceptions))
    (root / "liftoff.gff3").write_text(
        _gff("Liftoff", 1000, 1000, 1050, 2001, 2201, with_exceptions))
    # miniprot found gene A only after its selenocysteine: residues 4-26, the
    # fragment that beat the truncated score of the whole model (SEPHS2).
    (root / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t1110\t1328\t.\t+\t.\tID=MP1;Rank=1;Target=txA 4 26\n"
        "chr1\tminiprot\tCDS\t1110\t1160\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t1311\t1328\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n")
    return root


@pytest.fixture
def hermetic_pipeline(monkeypatch):
    from lifton import lifton_utils as _lu, run_liftoff, run_miniprot
    monkeypatch.setattr(_lu, "check_miniprot_installed", lambda: None)
    monkeypatch.setattr(run_miniprot, "check_miniprot_installed", lambda: True)

    def _fail(*args, **kwargs):
        raise AssertionError("External liftoff/miniprot must NOT be invoked in tests")

    monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)
    monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)


def _run(work, name, *extra):
    from lifton import lifton as lifton_main
    out = work / f"{name}.gff3"
    argv = [str(work / "tgt.fa"), str(work / "ref.fa"),
            "-g", str(work / "ref.gff3"), "-L", str(work / "liftoff.gff3"),
            "-M", str(work / "miniprot.gff3"), "-o", str(out),
            "-ad", "RefSeq", "--force", "-dir", str(work / f"dir_{name}"), *extra]
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    manifest = json.loads((work / f"dir_{name}" / "run_manifest.json").read_text())
    return out.read_text(), manifest


def _rows(text, parent=None, feature_id=None, ftype=None):
    out = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        attrs = dict(f.split("=", 1) for f in cols[8].split(";") if "=" in f)
        if ftype and cols[2] != ftype:
            continue
        if parent and attrs.get("Parent") != parent:
            continue
        if feature_id and attrs.get("ID") != feature_id:
            continue
        out.append((cols, attrs))
    return out


class TestEndToEnd:
    def test_the_selenoprotein_keeps_its_whole_lifted_model(self, tmp_path, hermetic_pipeline):
        text, _ = _run(_workspace(tmp_path / "w"), "out")
        (mrna, attrs), = _rows(text, feature_id="txA", ftype="mRNA")
        assert attrs["protein_identity"] == "1.000"
        assert "mutation" not in attrs
        cds = [(int(c[3]), int(c[4])) for c, _ in _rows(text, parent="txA", ftype="CDS")]
        assert cds == [(1101, 1160), (1311, 1328)], "the N-terminus was cut off"

    def test_transl_except_is_rewritten_into_target_coordinates(self, tmp_path, hermetic_pipeline):
        text, manifest = _run(_workspace(tmp_path / "w"), "out")
        a = {a.get("transl_except") for _, a in _rows(text, parent="txA", ftype="CDS")}
        b = {a.get("transl_except") for _, a in _rows(text, parent="txB", ftype="CDS")}
        assert a == {"(pos:1107..1109%2Caa:Sec)"}
        assert b == {"(pos:complement(join(1559..1560%2C1711))%2Caa:Sec)"}
        assert "107..109" not in text and "559..560" not in text, \
            "a reference coordinate reached the output"
        assert manifest["counts"]["transl_except_written"] == 2

    def test_a_target_that_codes_cysteine_carries_no_exception(self, tmp_path, hermetic_pipeline):
        text, manifest = _run(_workspace(tmp_path / "w"), "out")
        for _, attrs in _rows(text, parent="txC", ftype="CDS"):
            assert "transl_except" not in attrs
        assert manifest["counts"]["transl_except_not_applicable"] == 1

    def test_the_minus_strand_split_codon_scores_one(self, tmp_path, hermetic_pipeline):
        text, _ = _run(_workspace(tmp_path / "w"), "out")
        (_, attrs), = _rows(text, feature_id="txB", ftype="mRNA")
        assert attrs["protein_identity"] == "1.000"

    def test_a_gene_without_exceptions_is_byte_identical(self, tmp_path, hermetic_pipeline):
        declared, _ = _run(_workspace(tmp_path / "w"), "out")
        control, _ = _run(_workspace(tmp_path / "c", with_exceptions=False), "out")

        def gene_d(text):
            return [line for line in text.splitlines() if "D;" in line or "=txD" in line
                    or "ID=geneD" in line]
        assert gene_d(declared) == gene_d(control) != []

    def test_threads_agree(self, tmp_path, hermetic_pipeline):
        serial, _ = _run(_workspace(tmp_path / "w"), "serial")
        threaded, _ = _run(tmp_path / "w", "threaded", "-t", "2", "--locus-pipeline")
        assert serial == threaded

    def test_supplied_proteins_give_the_same_result(self, tmp_path, hermetic_pipeline):
        work = _workspace(tmp_path / "w")
        first, _ = _run(work, "extracted")
        intermediate = work / "dir_extracted"
        proteins = next(intermediate.rglob("proteins.fa"))
        transcripts = next(intermediate.rglob("transcripts.fa"))
        supplied, _ = _run(work, "supplied", "-P", str(proteins), "-T", str(transcripts))
        assert supplied == first

    def test_an_ncbi_protein_fasta_with_u_gives_the_same_result(self, tmp_path, hermetic_pipeline):
        work = _workspace(tmp_path / "w")
        first, _ = _run(work, "extracted")
        intermediate = work / "dir_extracted"
        proteins = next(intermediate.rglob("proteins.fa"))
        transcripts = next(intermediate.rglob("transcripts.fa"))
        # NCBI protein FASTAs write selenocysteine as U.
        records = proteins.read_text().split(">")[1:]
        ncbi = work / "ncbi_proteins.fa"
        with ncbi.open("w") as handle:
            for record in records:
                name, _, seq = record.partition("\n")
                seq = seq.replace("\n", "")
                body, end = (seq[:-1], "*") if seq.endswith("*") else (seq, "")
                handle.write(f">{name}\n{body.replace('*', 'U')}{end}\n")
        supplied, _ = _run(work, "ncbi", "-P", str(ncbi), "-T", str(transcripts))
        assert _rows(supplied, feature_id="txA", ftype="mRNA")[0][1]["protein_identity"] == "1.000"
        assert ({a.get("transl_except") for _, a in _rows(supplied, parent="txB", ftype="CDS")}
                == {"(pos:complement(join(1559..1560%2C1711))%2Caa:Sec)"})


# ---------------------------------------------------------------------------
# The two paths the end-to-end fixture cannot reach (its DNA is identical)
# ---------------------------------------------------------------------------

class TestVariantCall:
    def _trans(self):
        from lifton import lifton_class
        return lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)

    def test_a_declared_stop_is_not_a_gained_stop(self):
        from lifton import lifton_class, variants
        status = lifton_class.Lifton_Status()
        aln, peps = self._trans().align_coding_seq(
            "MA*AAK*", "MA*AAR*", status, readthrough=frozenset({2}))
        assert peps == ["MAUAAK", ""]
        dna = SimpleNamespace(identity=0.99, query_aln="ATGGCTTGAGCTGCTAAATAA",
                              ref_aln="ATGGCTTGAGCTGCTAGATAA")
        variants.find_variants(dna, aln, status, peps, False)
        assert "stop_codon_gain" not in status.status

    def test_without_a_declaration_it_still_is(self):
        from lifton import lifton_class
        status = lifton_class.Lifton_Status()
        _, peps = self._trans().align_coding_seq("MA*AAK*", "MA*AAR*", status)
        assert len(peps) == 3


class TestStopCompletion:
    def test_completing_onto_a_declared_codon_is_refused(self):
        from lifton import orf_completion
        ends_at_sec = SimpleNamespace(query_aln="MA*---", ref_aln="MA*AA*")
        whole = SimpleNamespace(query_aln="MAUAA*", ref_aln="MA*AA*")
        assert orf_completion._ends_on_readthrough(ends_at_sec, frozenset({2}))
        assert not orf_completion._ends_on_readthrough(whole, frozenset({2}))
        assert not orf_completion._ends_on_readthrough(ends_at_sec, None)


class TestEvaluator:
    def test_every_tool_is_scored_through_the_declared_codon(self, tmp_path, hermetic_pipeline):
        """The benchmark evaluator imports LiftOn's alignment, so it rated an
        identical selenoprotein as truncated -- for Liftoff, miniprot and
        LiftOn alike. It now reads the declaration from the REFERENCE."""
        import pyfaidx
        from benchmarks.compare import evaluator
        work = _workspace(tmp_path / "w")
        _run(work, "out")
        ref, _ = evaluator.build_reference(str(work / "ref.gff3"), str(work / "ref.fa"),
                                           log=lambda *args: None)
        assert ref["txA"]["readthrough"] == frozenset({2})
        assert ref["txD"]["readthrough"] == frozenset()
        db = evaluator._build_db(str(work / "out.gff3"))
        fasta = pyfaidx.Fasta(str(work / "tgt.fa"))
        for tx in ("txA", "txB"):
            mrna = db[tx]
            record = evaluator._eval_one_mrna(
                mrna, evaluator._children(db, mrna, "exon"),
                evaluator._children(db, mrna, ("CDS", "stop_codon")), tx, ref, fasta, False)
            assert record["protein_identity"] == 1.0, tx
            assert record["orf_valid"] == 1, tx
