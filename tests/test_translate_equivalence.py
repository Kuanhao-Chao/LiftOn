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
