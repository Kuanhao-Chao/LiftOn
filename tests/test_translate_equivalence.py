"""The codon-map translation must equal Biopython's, on every input.

Profiling Step 7 on whole genomes found `Bio.Data.CodonTable.__getitem__`
called 62 million times on rice and 257 million on dog to cat -- one
Python-level call per codon -- inside a translation block that is 8-9 % of the
phase (notes/step7_profile_2026-09.md). `coding.translate` answers from a flat
64-entry dict instead, and sends anything that dict cannot answer to Biopython
untouched.

That makes equivalence provable rather than arguable: the fast path's domain is
64 codons, so it can be checked exhaustively, and everything outside it is
Biopython by construction. These tests do both.
"""
from __future__ import annotations

import itertools
import random

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

from lifton import coding


BASES = "ACGT"
CODONS = ["".join(c) for c in itertools.product(BASES, repeat=3)]
TABLES = sorted(CodonTable.unambiguous_dna_by_id)


def _biopython(dna, table_id):
    complete = len(dna) - len(dna) % 3
    if not complete:
        return ""
    return str(Seq(dna[:complete]).translate(table=table_id))


class TestExhaustive:
    @pytest.mark.parametrize("table_id", TABLES)
    def test_every_codon_of_every_ncbi_table(self, table_id):
        # The fast path's entire domain, one table at a time.
        for codon in CODONS:
            assert coding.translate(codon, table_id) == \
                _biopython(codon, table_id), (codon, table_id)

    @pytest.mark.parametrize("table_id", TABLES)
    def test_the_map_covers_all_64_codons(self, table_id):
        codon_map = coding._codon_map(table_id)
        assert set(codon_map) == set(CODONS)
        table = CodonTable.unambiguous_dna_by_id[table_id]
        # A codon that is BOTH a stop and an amino acid (tables 27, 28, 31)
        # translates as the amino acid, which is what Biopython does. So the
        # "*" entries are the stops that carry no amino acid -- not every stop.
        assert {c for c, aa in codon_map.items() if aa == "*"} == \
            set(table.stop_codons) - set(table.forward_table)

    @pytest.mark.parametrize("table_id", [27, 28, 31])
    def test_a_dual_meaning_codon_translates_as_its_amino_acid(self, table_id):
        table = CodonTable.unambiguous_dna_by_id[table_id]
        dual = set(table.stop_codons) & set(table.forward_table)
        assert dual, f"table {table_id} was chosen for having dual codons"
        for codon in dual:
            assert coding.translate(codon, table_id) == \
                table.forward_table[codon] != "*"

    def test_the_orf_scan_still_treats_a_dual_codon_as_a_stop(self):
        # stop_codons answers a different question from translate: for finding
        # an open reading frame, a codon the table lists as a stop is a stop.
        assert "TGA" in coding.stop_codons(27)

    def test_a_whole_sequence_of_every_codon(self):
        sequence = "".join(CODONS)
        for table_id in TABLES:
            assert coding.translate(sequence, table_id) == \
                _biopython(sequence, table_id)


class TestOutsideTheFastPath:
    @pytest.mark.parametrize("sequence", [
        "ATGNNNTAA",            # N runs -> X
        "ATGRYKTAA",            # IUPAC ambiguity
        "atgtgaagataa",         # soft-masked reference sequence
        "AtGtGaAgAtAa",         # mixed case
        "ATG-TG-TAA",           # gap characters
        "ATGTGAAGATA",          # trailing partial codon
        "AT",                   # shorter than one codon
        "",
    ])
    @pytest.mark.parametrize("table_id", [1, 2, 11])
    def test_biopython_answers_and_we_agree(self, sequence, table_id):
        try:
            expected = _biopython(sequence, table_id)
        except Exception as error:                      # noqa: BLE001
            with pytest.raises(type(error)):
                coding.translate(sequence, table_id)
            return
        assert coding.translate(sequence, table_id) == expected

    def test_lowercase_is_not_quietly_mistranslated(self):
        # The reference genome is soft-masked, so get_dna_sequence hands us
        # mixed case. Getting this wrong would silently corrupt those proteins.
        for table_id in (1, 2, 11):
            assert coding.translate("atgtgaagataa", table_id) == \
                coding.translate("ATGTGAAGATAA", table_id)


class TestFuzz:
    @pytest.mark.parametrize("table_id", [1, 2, 4, 5, 11, 33])
    def test_random_sequences_agree(self, table_id):
        rng = random.Random(20260919 + table_id)
        alphabet = BASES + BASES.lower() + "N"
        for _ in range(300):
            length = rng.randrange(0, 240)
            sequence = "".join(rng.choice(alphabet) for _ in range(length))
            assert coding.translate(sequence, table_id) == \
                _biopython(sequence, table_id), sequence


class TestEscapeHatch:
    def test_legacy_translate_reproduces_the_same_answer(self, monkeypatch):
        sequence = "".join(CODONS)
        monkeypatch.setenv("LIFTON_LEGACY_TRANSLATE", "1")
        legacy = coding.translate(sequence, 2)
        monkeypatch.delenv("LIFTON_LEGACY_TRANSLATE")
        assert legacy == coding.translate(sequence, 2) == _biopython(sequence, 2)


class TestAttributeEncoding:
    """The GFF3 attribute encoder answers from `str.translate` behind a guard
    that returns the string untouched when nothing needs encoding.

    5.8 M calls on rice and 17.3 M on dog to cat (5 % and 3 % of Step 7) ran a
    Python-level loop over every character of every value, almost always to
    hand it back unchanged.
    """

    @staticmethod
    def _legacy(value):
        """The per-character loop this replaced, kept verbatim as the oracle."""
        from lifton.io.ncbi_gff3_spec import RESERVED_CHARS
        if value is None:
            return ""
        text = str(value)
        if not text:
            return text
        out = []
        for ch in text:
            if ch == "%":
                out.append("%25")
            elif ch in RESERVED_CHARS:
                out.append("%{:02X}".format(ord(ch)))
            else:
                out.append(ch)
        return "".join(out)

    def test_every_code_point_a_value_can_hold(self):
        from lifton.io.gff3_writer import encode_attribute_value
        for code_point in range(0x3000):
            character = chr(code_point)
            assert encode_attribute_value(character) == \
                self._legacy(character), hex(code_point)

    @pytest.mark.parametrize("value", [
        None, "", "GeneID:1,HGNC:HGNC:2", "a%00b", "%25", "tab\there",
        "semi;colon", "eq=uals", "amp&ersand", "new\nline", "cr\rlf",
        "%already%encoded%", "no-reserved-characters-at-all",
    ])
    def test_edge_cases(self, value):
        from lifton.io.gff3_writer import encode_attribute_value
        assert encode_attribute_value(value) == self._legacy(value)

    def test_fuzz(self):
        from lifton.io.gff3_writer import encode_attribute_value
        rng = random.Random(20260919)
        alphabet = ("ABCXYZabcxyz0189_.:-|%;=&,\t\n\r")
        for _ in range(5000):
            value = "".join(rng.choice(alphabet)
                            for _ in range(rng.randrange(0, 60)))
            assert encode_attribute_value(value) == self._legacy(value), value

    def test_a_clean_value_is_returned_unchanged(self):
        from lifton.io.gff3_writer import encode_attribute_value
        value = "rna-NM_001303012.2"
        assert encode_attribute_value(value) is value


class TestCloneAttributes:
    """Cloning a gffutils Attributes copies its backing mapping directly.

    ``items()`` builds a list and calls ``__getitem__`` per key; ``__setitem__``
    then re-wraps each value. Profiling Step 7 on dog to cat counted 85.5
    million ``__setitem__`` calls behind 8.85 million clones.
    """

    @staticmethod
    def _legacy(attributes):
        if attributes is None:
            return None
        clone = attributes.__class__()
        for key, value in attributes.items():
            clone[key] = list(value) if isinstance(value, (list, tuple)) else value
        return clone

    def _attributes(self, pairs):
        from gffutils.attributes import Attributes
        attributes = Attributes()
        for key, value in pairs:
            attributes[key] = value
        return attributes

    @pytest.mark.parametrize("pairs", [
        [],
        [("ID", ["cds-NP_001.1"])],
        [("ID", ["c"]), ("Parent", ["r"]),
         ("Dbxref", ["GeneID:1", "HGNC:HGNC:2"]), ("product", ["a protein"])],
        [("single", "not-a-list")],
        [("tuple", ("a", "b"))],
    ])
    def test_it_matches_the_loop_it_replaced(self, pairs):
        from lifton import coreutils
        attributes = self._attributes(pairs)
        assert dict(coreutils.clone_attributes(attributes)._d) == \
            dict(self._legacy(attributes)._d)

    def test_the_clone_is_independent_of_its_source(self):
        from lifton import coreutils
        attributes = self._attributes([("Dbxref", ["GeneID:1", "HGNC:2"])])
        clone = coreutils.clone_attributes(attributes)
        clone["Dbxref"].append("added")
        clone["ID"] = ["new"]
        assert attributes["Dbxref"] == ["GeneID:1", "HGNC:2"]
        assert "ID" not in attributes

    def test_the_class_is_preserved(self):
        from gffutils.attributes import Attributes
        from lifton import coreutils
        assert isinstance(coreutils.clone_attributes(self._attributes([])),
                          Attributes)
        assert type(coreutils.clone_attributes({"a": ["b"]})) is dict

    def test_none_and_plain_mappings_are_unchanged(self):
        from lifton import coreutils
        assert coreutils.clone_attributes(None) is None
        plain = {"ID": ["x"], "count": 1}
        clone = coreutils.clone_attributes(plain)
        assert clone == plain and clone is not plain
        clone["ID"].append("y")
        assert plain["ID"] == ["x"]
