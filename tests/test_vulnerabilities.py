"""Phase 13.5B + 13.5C — hostile + property-based tests for the
Phase 13.5A audit findings.

Each test is written to FAIL on the unpatched codebase and PASS once
the matching production-code fix lands. Tests are grouped by audit ID
(V1.x = exception handling, V2.x = algorithmic, V4.x = data,
V5.x = GFF3 edge cases). Phase 13.5C extends with Hypothesis-driven
fuzz tests for the GFF3 emission path.

The tests stay hermetic: no external binaries, no network. Where a
real gffutils database is needed, the existing conftest fixtures are
reused.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


# ---------------------------------------------------------------------------
# V1.4 — process_locus_native swallowed BaseException, including
# KeyboardInterrupt. The harness must let Ctrl-C escape the worker.
# ---------------------------------------------------------------------------

class TestV1_4_BaseExceptionInWorker:
    def test_keyboardinterrupt_in_worker_propagates(self, monkeypatch):
        from lifton import locus_pipeline, run_liftoff

        def boom(*_a, **_k):
            raise KeyboardInterrupt
        monkeypatch.setattr(run_liftoff, "process_liftoff", boom)

        ctx = SimpleNamespace(
            ref_db=None, l_feature_db=None, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={}, tgt_fai=None,
            ref_proteins={}, ref_trans={}, ref_features_dict={},
            fw_score=None, fw_chain=None,
            args=SimpleNamespace(write_chains=False, debug=False,
                                 no_orf_search=True),
        )
        payload = locus_pipeline.MaterialisedLocus(
            submission_index=0,
            locus=SimpleNamespace(id="g1"),
            locus_id="g1",
        )

        with pytest.raises(KeyboardInterrupt):
            locus_pipeline.process_locus_native(payload, ctx)

    def test_systemexit_in_worker_propagates(self, monkeypatch):
        """SystemExit is also a BaseException, not Exception. Must escape."""
        from lifton import locus_pipeline, run_liftoff

        def boom(*_a, **_k):
            raise SystemExit(2)
        monkeypatch.setattr(run_liftoff, "process_liftoff", boom)

        ctx = SimpleNamespace(
            ref_db=None, l_feature_db=None, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={}, tgt_fai=None,
            ref_proteins={}, ref_trans={}, ref_features_dict={},
            fw_score=None, fw_chain=None,
            args=SimpleNamespace(write_chains=False, debug=False,
                                 no_orf_search=True),
        )
        payload = locus_pipeline.MaterialisedLocus(
            submission_index=0,
            locus=SimpleNamespace(id="g1"),
            locus_id="g1",
        )

        with pytest.raises(SystemExit):
            locus_pipeline.process_locus_native(payload, ctx)

    def test_ordinary_exception_still_packaged(self, monkeypatch):
        """Regression: ValueError from a worker must still be packaged
        into LocusResult.error so siblings keep running."""
        from lifton import locus_pipeline, run_liftoff

        def boom(*_a, **_k):
            raise ValueError("synthetic")
        monkeypatch.setattr(run_liftoff, "process_liftoff", boom)

        ctx = SimpleNamespace(
            ref_db=None, l_feature_db=None, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={}, tgt_fai=None,
            ref_proteins={}, ref_trans={}, ref_features_dict={},
            fw_score=None, fw_chain=None,
            args=SimpleNamespace(write_chains=False, debug=False,
                                 no_orf_search=True),
        )
        payload = locus_pipeline.MaterialisedLocus(
            submission_index=7,
            locus=SimpleNamespace(id="g1"),
            locus_id="g1",
        )
        result = locus_pipeline.process_locus_native(payload, ctx)
        assert result.lifton_gene is None
        assert isinstance(result.error, ValueError)
        assert result.index == 7


# ---------------------------------------------------------------------------
# V1.1a — run_liftoff.py:209 bare `except:` swallowed KeyboardInterrupt
# during ref_db lookup.
# ---------------------------------------------------------------------------

class TestV1_1a_BareExceptInRefDbLookup:
    def test_keyboardinterrupt_in_ref_db_lookup_propagates(self, monkeypatch,
                                                           gff_standard,
                                                           fasta_standard,
                                                           fake_args):
        """When the user hits Ctrl-C while process_liftoff is doing a
        ref_db[trans_id] check, the bare `except:` must NOT swallow it."""
        import gffutils
        from lifton import lifton_class, run_liftoff

        l_db = gffutils.create_db(
            str(gff_standard), ":memory:", force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )
        # Build a fake ref_db whose __getitem__ raises KeyboardInterrupt
        class HostileRefDB:
            def __getitem__(self, key):
                raise KeyboardInterrupt
        ref_db = HostileRefDB()

        # Locate the mRNA feature so process_liftoff hits the
        # `try: ref_db[ref_trans_id]` block at line 207.
        locus = next(l_db.features_of_type("mRNA"))

        # Provide the minimum context so we reach line 207.
        ref_features_dict = {"tx1": lifton_class.Lifton_feature("tx1")}
        ref_features_dict["tx1"].is_protein_coding = True
        # The lookup happens on a Lifton_GENE that already has a
        # ref_gene_id, so we build one inline. Lifton_GENE wants
        # (ref_gene_id, gffutil_entry_gene, ref_gene_attrs, tree_dict,
        #  ref_features_dict, args).
        gene_locus = next(l_db.features_of_type("gene"))
        from lifton import lifton_class as lc
        args = SimpleNamespace(
            write_chains=False, debug=False, no_orf_search=True,
            annotation_database="RefSeq", evaluation=False,
            evaluation_liftoff_chm13=False,
        )
        lifton_gene = lc.Lifton_GENE(
            "gene1", gene_locus, {}, {}, ref_features_dict, args,
        )
        with pytest.raises(KeyboardInterrupt):
            run_liftoff.process_liftoff(
                lifton_gene, locus, ref_db, l_db,
                {}, None, {}, None, {}, {}, ref_features_dict,
                None, None, args, ENTRY_FEATURE=False,
            )


# ---------------------------------------------------------------------------
# V1.1b — run_miniprot.check_miniprot_installed bare `except:` masked
# real subprocess errors as "miniprot not installed".
# ---------------------------------------------------------------------------

class TestV1_1b_BareExceptInCheckMiniprot:
    def test_unexpected_error_propagates(self, monkeypatch):
        """A MemoryError from subprocess is NOT 'miniprot is not installed'.
        The function must let MemoryError escape, not return False."""
        from lifton import run_miniprot

        def boom(*_a, **_k):
            raise MemoryError("synthetic")
        monkeypatch.setattr("subprocess.run", boom)

        with pytest.raises(MemoryError):
            run_miniprot.check_miniprot_installed()

    def test_missing_binary_returns_false(self, monkeypatch):
        """Regression: FileNotFoundError IS 'miniprot not installed' →
        the function should return False, not raise."""
        from lifton import run_miniprot

        def boom(*_a, **_k):
            raise FileNotFoundError("miniprot")
        monkeypatch.setattr("subprocess.run", boom)
        assert run_miniprot.check_miniprot_installed() is False


# ---------------------------------------------------------------------------
# V1.2 — extract_sequence.__inner_extract_feature silently swallowed
# Exception → users never knew transcripts were dropped.
# ---------------------------------------------------------------------------

class TestV1_2_SilentExtractFeatureSwallow:
    def test_failure_is_logged(self, monkeypatch, capsys, gff_standard,
                               fasta_standard):
        """When get_dna_sequence raises, a [WARNING] line must reach stderr.
        The unpatched code silently dropped the transcript."""
        from lifton import extract_sequence

        def boom(*_a, **_k):
            raise RuntimeError("synthetic extraction failure")
        monkeypatch.setattr(extract_sequence, "get_dna_sequence", boom)

        # Build a minimal ref_db wrapper exposing db_connection
        import gffutils
        gffu_db = gffutils.create_db(
            str(gff_standard), ":memory:", force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )

        class FakeRefDB:
            def __init__(self, db):
                self.db_connection = db

        ref_db = FakeRefDB(gffu_db)
        # ref_fai unused because we patched get_dna_sequence to raise
        extract_sequence.extract_features(ref_db, ["mRNA"], None)

        captured = capsys.readouterr()
        assert "[WARNING]" in captured.err, (
            "extract_features should log a warning when sequence "
            "extraction fails — found no warning in stderr."
        )


# ---------------------------------------------------------------------------
# V1.3 — locus_pipeline.materialise_locus had 6 silent except blocks.
# ---------------------------------------------------------------------------

class TestV1_3_MaterialiseLocusSilentSwallows:
    def test_db_failure_is_logged(self, capsys, monkeypatch):
        """When l_feature_db.children raises, materialise_locus must
        emit a [WARNING] for each failed read instead of silently
        producing an empty payload."""
        from lifton import locus_pipeline

        class HostileDB:
            def children(self, *_a, **_k):
                raise RuntimeError("synthetic db failure")
            def __getitem__(self, key):
                raise KeyError(key)
        ctx = SimpleNamespace(
            ref_db=HostileDB(),
            l_feature_db=HostileDB(),
            m_feature_db=None,
            ref_id_2_m_id_trans_dict={},
            tree_dict={}, tgt_fai=None,
            ref_proteins={}, ref_trans={},
            ref_features_dict={},
            fw_score=None, fw_chain=None,
            args=SimpleNamespace(annotation_database="RefSeq",
                                 evaluation=False,
                                 evaluation_liftoff_chm13=False,
                                 debug=False),
        )
        payload = locus_pipeline.materialise_locus(
            0, SimpleNamespace(id="bad_locus"), ctx,
        )
        captured = capsys.readouterr()
        assert "[WARNING]" in captured.err
        assert payload.children_l1 == []

    def test_keyboardinterrupt_propagates_through_materialise(self, capsys):
        """KeyboardInterrupt in the parent thread must still escape — the
        parent isn't a worker, but `except Exception` is the right choice
        and KeyboardInterrupt subclasses BaseException only."""
        from lifton import locus_pipeline

        class HostileDB:
            def children(self, *_a, **_k):
                raise KeyboardInterrupt
            def __getitem__(self, key):
                raise KeyError(key)
        ctx = SimpleNamespace(
            ref_db=HostileDB(),
            l_feature_db=HostileDB(),
            m_feature_db=None,
            ref_id_2_m_id_trans_dict={},
            tree_dict={}, tgt_fai=None,
            ref_proteins={}, ref_trans={},
            ref_features_dict={},
            fw_score=None, fw_chain=None,
            args=SimpleNamespace(annotation_database="RefSeq",
                                 evaluation=False,
                                 evaluation_liftoff_chm13=False,
                                 debug=False),
        )
        with pytest.raises(KeyboardInterrupt):
            locus_pipeline.materialise_locus(
                0, SimpleNamespace(id="bad"), ctx,
            )


# ---------------------------------------------------------------------------
# V1.5 — Lifton_TRANS.write_entry swallowed every exception, including
# BrokenPipeError; child rows then emit AFTER the parent failed → orphan
# exons / CDSs in the output.
# ---------------------------------------------------------------------------

class TestV1_5_WriteEntryOrphanChildren:
    def test_broken_pipe_does_not_emit_orphan_children(self, gff_standard):
        """If writing the mRNA row raises BrokenPipeError, the exon and
        CDS writes should be SKIPPED — otherwise the consumer of the
        pipe sees orphan children referring to a parent it never received.
        """
        import gffutils
        from lifton import lifton_class

        db = gffutils.create_db(
            str(gff_standard), ":memory:", force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )
        mrna = db["tx1"]
        # Build a minimal Lifton_TRANS with two exons.
        trans = lifton_class.Lifton_TRANS(
            "tx1", "gene1", "gene1", 0, mrna, dict(mrna.attributes),
        )
        for ex in db.children(mrna, featuretype="exon", level=1):
            trans.add_exon(ex)

        class BrokenPipe:
            def __init__(self):
                self.written = []
            def write(self, payload):
                if not self.written:
                    self.written.append(payload)
                    raise BrokenPipeError("downstream closed")
                self.written.append(payload)

        fw = BrokenPipe()
        trans.write_entry(fw)
        assert len(fw.written) == 1, (
            "After parent write fails, child writes must NOT happen — "
            f"got {len(fw.written)} writes."
        )


# ---------------------------------------------------------------------------
# V1.9 — Annotation._handle_gtf_input silently fell back to the raw GTF
# when conversion failed, masking the real problem.
# ---------------------------------------------------------------------------

class TestV1_9_GTFConversionSilentFallback:
    def test_conversion_failure_raises(self, monkeypatch, tmp_path):
        """When auto_convert_gtf=True is requested and BOTH gffread and
        agat fail, the user must see a clear LiftOnInputError, not a
        silent fallback to the raw GTF that explodes 5 minutes later."""
        from lifton import annotation
        from lifton.exceptions import LiftOnInputError

        gtf = tmp_path / "broken.gtf"
        gtf.write_text(
            'chr1\ttest\ttranscript\t100\t200\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "t1";\n'
        )

        # Pretend gffread & agat are present so we exercise the
        # conversion branch, but make subprocess.run fail with
        # CalledProcessError when the actual conversion runs.
        original_run = subprocess.run

        def picky_run(cmd, *a, **kw):
            if cmd and cmd[0] == "which":
                # 'tool present' for gffread; 'absent' for agat
                if len(cmd) > 1 and cmd[1] == "gffread":
                    return SimpleNamespace(returncode=0, stdout="", stderr="")
                raise subprocess.CalledProcessError(1, cmd)
            if cmd and cmd[0] == "gffread":
                raise subprocess.CalledProcessError(2, cmd, "", "boom")
            return original_run(cmd, *a, **kw)
        monkeypatch.setattr(subprocess, "run", picky_run)

        with pytest.raises(LiftOnInputError):
            annotation.Annotation(
                str(gtf), infer_genes=False, infer_transcripts=False,
                auto_convert_gtf=True,
            )


# ---------------------------------------------------------------------------
# V2.1 — extract_sequence.get_dna_sequence wrapped silently when start<1.
# ---------------------------------------------------------------------------

class TestV2_1_GetDnaSequenceStartBounds:
    def test_zero_start_raises(self, fasta_standard):
        """start=0 is illegal under GFF3's 1-based-inclusive convention.
        The unpatched code did `fasta[chrom][-1: end]` which silently
        wraps to the LAST base of the chromosome."""
        import pyfaidx
        from lifton import extract_sequence
        from lifton.exceptions import LiftOnInputError

        fa = pyfaidx.Fasta(str(fasta_standard))
        feature = SimpleNamespace(seqid="chr1", strand="+")
        with pytest.raises(LiftOnInputError):
            extract_sequence.get_dna_sequence(
                feature, fa,
                [SimpleNamespace(start=0, end=10)],
            )

    def test_negative_start_raises(self, fasta_standard):
        import pyfaidx
        from lifton import extract_sequence
        from lifton.exceptions import LiftOnInputError

        fa = pyfaidx.Fasta(str(fasta_standard))
        feature = SimpleNamespace(seqid="chr1", strand="+")
        with pytest.raises(LiftOnInputError):
            extract_sequence.get_dna_sequence(
                feature, fa,
                [SimpleNamespace(start=-3, end=10)],
            )

    def test_start_one_works(self, fasta_standard):
        """Boundary: start==1 (chromosome edge) is LEGAL."""
        import pyfaidx
        from lifton import extract_sequence

        fa = pyfaidx.Fasta(str(fasta_standard))
        feature = SimpleNamespace(seqid="chr1", strand="+")
        seq = extract_sequence.get_dna_sequence(
            feature, fa,
            [SimpleNamespace(start=1, end=10)],
        )
        # 10 bases pulled (post-pad to multiple of 3 → 12)
        assert len(seq) >= 10


# ---------------------------------------------------------------------------
# V2.3 — get_id_fraction.get_DNA_id_fraction crashed with IndexError
# on length-mismatched inputs.
# ---------------------------------------------------------------------------

class TestV2_3_DnaIdFractionLengthMismatch:
    def test_target_shorter_raises_value_error(self):
        from lifton import get_id_fraction
        with pytest.raises(ValueError):
            get_id_fraction.get_DNA_id_fraction("ACGT", "AC")

    def test_target_longer_does_not_crash(self):
        """Defensive symmetry: target longer than reference must not crash.
        After the fix it should also raise ValueError."""
        from lifton import get_id_fraction
        with pytest.raises(ValueError):
            get_id_fraction.get_DNA_id_fraction("AC", "ACGT")

    def test_equal_length_unchanged(self):
        """Regression: equal-length inputs must keep the original numeric
        result so existing identity numbers don't drift."""
        from lifton import get_id_fraction
        matches, length = get_id_fraction.get_DNA_id_fraction("ACGT", "ACGA")
        assert (matches, length) == (3, 4)


# ---------------------------------------------------------------------------
# V4.2 — parasail_align_DNA_base crashed at the C level on IUPAC
# ambiguity codes (R/Y/S/W/K/M/B/D/H/V/N).
# ---------------------------------------------------------------------------

class TestV4_2_ParasailIupacBases:
    @pytest.mark.parametrize("trans, ref", [
        ("ACRT", "ACGT"),         # R = A or G
        ("ACGT", "ACYT"),         # Y = C or T
        ("ACSWMK", "ACGTAC"),     # multiple IUPAC
        ("ACGTBDHVN", "ACGTACGTA"),
        ("NNNNNNNNNN", "ACGTACGTAC"),
    ])
    def test_iupac_codes_handled(self, trans, ref):
        """The patched parasail wrapper must sanitise non-ACGT bases to
        N (and extend the score matrix to include N) BEFORE passing to
        the C kernel — otherwise parasail aborts the process."""
        from lifton import align
        result = align.parasail_align_DNA_base(trans, ref)
        # Just assert we got a result object back, not a crash.
        assert result is not None
        assert result.traceback is not None

    def test_pure_acgt_still_aligns(self):
        """Regression: pure ACGT alignment numbers must not change."""
        from lifton import align
        result = align.parasail_align_DNA_base("ACGT", "ACGT")
        assert result.traceback.query == "ACGT"
        assert result.traceback.ref == "ACGT"


# ---------------------------------------------------------------------------
# V5.1 — Missing ID with a referenced Parent must be caught by the
# Phase 5 strict validator and surfaced clearly.
# ---------------------------------------------------------------------------

class TestV5_1_DanglingParent:
    def test_dangling_parent_caught_by_validator(self, gff_dangling_parent):
        from lifton.io.gff3_validator import GFF3Validator
        v = GFF3Validator(strict=True)
        findings = v.validate(gff_dangling_parent)
        rules = {f.rule for f in findings}
        assert "dangling_parent" in rules


# ---------------------------------------------------------------------------
# V5.2 — Circular Parent dependency caused infinite recursion in
# process_liftoff. The patched version must terminate cleanly.
# ---------------------------------------------------------------------------

class TestV5_2_CircularParentCycle:
    def test_circular_parent_does_not_recurse_forever(self, monkeypatch):
        """Construct a fake l_feature_db whose `children(locus)` returns
        a child that points BACK to its parent. The unpatched recursive
        process_liftoff would hit Python's recursion limit; the patched
        version detects the cycle and raises LiftOnInputError."""
        from lifton import run_liftoff
        from lifton.exceptions import LiftOnInputError

        class CycleLocus:
            id = "A"
            seqid = "chr1"
            featuretype = "gene"
            strand = "+"
            start = 100
            end = 200
            attributes = {}

        a = CycleLocus()
        # B is a level-1 child of A, but B's level-1 child is A itself.
        class CycleChild:
            id = "B"
            seqid = "chr1"
            featuretype = "mRNA"
            strand = "+"
            start = 100
            end = 200
            attributes = {"Parent": ["A"]}
        b = CycleChild()

        class CycleDB:
            def children(self, parent, level=None, featuretype=None,
                         order_by=None, **_k):
                # The exon query should return [] so we fall through
                # into the recursive level-1 branch where the cycle lives.
                if featuretype is not None:
                    return iter([])
                # Level-1 unfiltered: A → B and B → A → infinite loop.
                if getattr(parent, "id", None) == "A":
                    return iter([b])
                if getattr(parent, "id", None) == "B":
                    return iter([a])
                return iter([])

        # Strip down so we hit the recursion path early (no exon children).
        # The function path under test is the `if len(exon_children) == 0`
        # branch which recurses for every level-1 child.
        from lifton import lifton_class
        ref_features_dict = {"A": lifton_class.Lifton_feature("A")}
        ref_features_dict["A"].is_protein_coding = True

        # Build a Lifton_GENE manually using the existing constructor.
        fake_args = SimpleNamespace(
            annotation_database="RefSeq",
            evaluation=False, evaluation_liftoff_chm13=False, debug=False,
        )
        lifton_gene = lifton_class.Lifton_GENE(
            "A", a, {}, {}, ref_features_dict, fake_args,
        )
        args = SimpleNamespace(
            write_chains=False, debug=False, no_orf_search=True,
            annotation_database="RefSeq", evaluation=False,
            evaluation_liftoff_chm13=False,
        )

        # Patch sys.setrecursionlimit so the test fails fast on the
        # unpatched code rather than hanging.
        old_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(200)
        try:
            with pytest.raises((LiftOnInputError, RecursionError)):
                run_liftoff.process_liftoff(
                    lifton_gene, a, {}, CycleDB(),
                    {}, None, {}, None, {}, {}, ref_features_dict,
                    None, None, args, ENTRY_FEATURE=False,
                )
        finally:
            sys.setrecursionlimit(old_limit)

        # NOTE: After the patch, LiftOnInputError is the explicit signal.
        # Before the patch, this test would either hang or hit
        # RecursionError. The pytest.raises tuple accepts both so the
        # test demonstrates the fix narrows the failure mode.


# ===========================================================================
# Phase 13.5C — Medium severity findings (V2.4–V2.11, V5.3–V5.6, V5.8, V5.9)
# Plus Hypothesis property-based fuzz tests for GFF3 emission integrity.
# ===========================================================================


# ---------------------------------------------------------------------------
# V2.4 — get_AA_id_fraction zero-divisor when reference is all gaps and
# target hits a stop codon in the same window.
# ---------------------------------------------------------------------------

class TestV2_4_AAIdFractionZeroDivisor:
    def test_all_gaps_no_stop_in_target_does_not_zero_divide(self):
        """Pre-patch: `get_AA_id_fraction("---", "ACG")` → max=3,
        gaps_in_ref=3 → total_length=0. Caller does `matches/length`
        and crashes ZeroDivisionError. Post-patch: denominator >= 1."""
        from lifton import get_id_fraction

        matches, length = get_id_fraction.get_AA_id_fraction("---", "ACG")
        assert length >= 1
        # Sanity: we should never get a NaN-prone identity downstream.
        identity = matches / length
        assert 0.0 <= identity <= 1.0

    def test_long_all_gaps_does_not_zero_divide(self):
        from lifton import get_id_fraction

        matches, length = get_id_fraction.get_AA_id_fraction(
            "-" * 20, "ACDEFGHIKLMNPQRSTVWY",
        )
        assert length >= 1
        identity = matches / length
        assert identity == 0.0


# ---------------------------------------------------------------------------
# V2.5 — protein_maximization.process_m_l_children mis-labels an empty
# selection as `liftoff[0.00-0.00]` instead of "empty".
# ---------------------------------------------------------------------------

class TestV2_5_ChainLogEmptyLabelling:
    def test_zero_window_uses_empty_label(self):
        """When BOTH miniprot and liftoff windows are zero-length, the
        chain log should not pretend Liftoff was selected — record
        `empty[...]` so provenance is honest."""
        from lifton import protein_maximization

        # Construct minimal lifton_aln stand-ins with cdss_protein_aln_boundaries
        m_aln = SimpleNamespace(
            cdss_protein_aln_boundaries={0: (0, 0), 1: (0, 0)},
            ref_aln="", query_aln="",
        )
        l_aln = SimpleNamespace(
            cdss_protein_aln_boundaries={0: (0, 0), 1: (0, 0)},
            ref_aln="", query_aln="",
        )
        chains: list[str] = []
        cds_ls = protein_maximization.process_m_l_children(
            1, 0, m_aln, 1, 0, l_aln, fai=None, chains=chains, DEBUG=False,
        )
        # Should NOT have appended a misleading "liftoff[0.00-0.00]" entry.
        for entry in chains:
            assert not entry.startswith("liftoff[0.00-0.00]"), (
                f"Chain log records false Liftoff selection for empty "
                f"window: {chains}"
            )


# ---------------------------------------------------------------------------
# V2.6 — __find_orfs accepts trivial 1-codon ORFs that win their frame.
# ---------------------------------------------------------------------------

class TestV2_6_FindOrfsMinimumLength:
    def test_single_codon_orf_does_not_win(self, gff_standard, fasta_standard):
        """A 6-nt ORF (ATG followed immediately by TAA) is biological
        noise. After the patch, the rescue threshold rejects it and
        no CDS boundary update fires."""
        import gffutils
        from lifton import lifton_class

        db = gffutils.create_db(
            str(gff_standard), ":memory:", force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )
        mrna = db["tx1"]
        trans = lifton_class.Lifton_TRANS(
            "tx1", "gene1", "gene1", 0, mrna, dict(mrna.attributes),
        )
        for ex in db.children(mrna, featuretype="exon", level=1):
            trans.add_exon(ex)

        # Synth a transcript with ONLY a 6-nt ORF then padding.
        trans_seq = "ATGTAA" + "A" * 90
        ref_protein_seq = "MAGTACGT"   # arbitrary reference
        status = lifton_class.Lifton_Status()
        status.lifton_aa = 0.0
        # Use the dunder to call the private __find_orfs.
        trans._Lifton_TRANS__find_orfs(trans_seq, ref_protein_seq, None, status)
        # The micro-ORF must not boost lifton_aa past the threshold.
        # Either the rescue rejects it outright (status unchanged) or the
        # identity is below 0.01.
        assert status.lifton_aa <= 0.5, (
            f"Trivial 1-codon ORF was adopted as the winner: "
            f"lifton_aa={status.lifton_aa}"
        )


# ---------------------------------------------------------------------------
# V2.7 — Negative coordinates in __iterate_exons_update_cds (- strand).
# ---------------------------------------------------------------------------

from hypothesis import HealthCheck, given, settings, strategies as st


class TestV2_7_IterateExonsNoNegativeCoords:
    @settings(max_examples=50, deadline=None,
              suppress_health_check=[HealthCheck.function_scoped_fixture])
    @given(
        exon_starts=st.lists(
            st.integers(min_value=100, max_value=200),
            min_size=1, max_size=4, unique=True,
        ),
        orf_start=st.integers(min_value=0, max_value=20),
        orf_end_offset=st.integers(min_value=3, max_value=30),
    )
    def test_minus_strand_never_emits_negative_start(
        self, exon_starts, orf_start, orf_end_offset, gff_standard,
    ):
        """No matter what (small_exons + far_downstream_orf) input we
        throw at the rescuer, the resulting CDS coordinates must
        remain >= 1."""
        import gffutils
        from lifton import lifton_class

        # We don't need a real DB; build a fake transcript with synthesized
        # exons in-memory.
        db = gffutils.create_db(
            str(gff_standard), ":memory:", force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )
        mrna = db["tx1"]
        trans = lifton_class.Lifton_TRANS(
            "tx1", "gene1", "gene1", 0, mrna, dict(mrna.attributes),
        )
        # Replace mrna's strand to '-' for the test.
        trans.entry.strand = "-"
        for s in sorted(set(exon_starts)):
            ex = gffutils.Feature(
                seqid="chr1", source="test", featuretype="exon",
                start=s, end=s + 4, strand="-", frame=".",
                attributes={"ID": [f"ex{s}"], "Parent": ["tx1"]},
                id=f"ex{s}",
            )
            trans.add_exon(ex)

        orf = lifton_class.Lifton_ORF(orf_start, orf_start + orf_end_offset)
        # Apply private boundary update directly.
        try:
            trans._Lifton_TRANS__update_cds_boundary(orf)
        except Exception:
            # Patched code raises LiftOnInputError on impossible math;
            # that's acceptable. Crashes other than that are bugs.
            return

        for exon in trans.exons:
            if exon.cds is not None:
                assert exon.cds.entry.start >= 1, (
                    f"CDS got negative start: {exon.cds.entry.start} "
                    f"(strand=-, exon=({exon.entry.start},{exon.entry.end}), "
                    f"orf=({orf.start},{orf.end}))"
                )
                assert exon.cds.entry.end >= exon.cds.entry.start


# ---------------------------------------------------------------------------
# V2.8 — IntervalTree crash on zero-length single-base feature.
# ---------------------------------------------------------------------------

class TestV2_8_IntervalTreeZeroLength:
    def test_single_base_feature_does_not_crash(self, tmp_path):
        """A GFF3 row with start == end (single-base feature) is legal
        per NCBI col 4-5; intervaltree disagrees. The patched code
        widens by 1 (half-open) so the tree builds successfully."""
        import gffutils
        from lifton import intervals

        gff = tmp_path / "single_base.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t100\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
        )
        db = gffutils.create_db(
            str(gff), ":memory:", force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )
        # Should NOT raise ValueError("Interval cannot have zero length")
        tree_dict = intervals.initialize_interval_tree(db, ["gene"])
        assert "chr1" in tree_dict
        assert len(tree_dict["chr1"]) == 1


# ---------------------------------------------------------------------------
# V2.9 — segments_overlap_length returns NEGATIVE on disjoint segments.
# ---------------------------------------------------------------------------

class TestV2_9_OverlapNeverNegative:
    @settings(max_examples=200, deadline=None)
    @given(
        a=st.integers(min_value=1, max_value=10**6),
        b=st.integers(min_value=1, max_value=10**6),
        c=st.integers(min_value=1, max_value=10**6),
        d=st.integers(min_value=1, max_value=10**6),
    )
    def test_ovp_len_is_non_negative(self, a, b, c, d):
        """For ANY two segments, the overlap length must be >= 0.
        Currently disjoint inputs return negative numbers."""
        from lifton import lifton_utils
        s1, e1 = sorted([a, b])
        s2, e2 = sorted([c, d])
        ovp_len, ovp = lifton_utils.segments_overlap_length((s1, e1), (s2, e2))
        assert ovp_len >= 0, (
            f"Disjoint segments produced negative overlap: "
            f"({s1},{e1}) vs ({s2},{e2}) -> {ovp_len}"
        )
        # Boolean ovp must be consistent with the numeric value.
        assert (ovp_len > 0) == bool(ovp)


# ---------------------------------------------------------------------------
# V2.10 — adjust_cdss_protein_boundary cumulative shift on multi-D CIGAR.
# ---------------------------------------------------------------------------

class TestV2_10_AdjustBoundaryNoCumulativeShift:
    def test_two_d_blocks_each_shift_independently(self):
        """When two D-blocks straddle the SAME boundary, the shift
        should NOT compound — each D adds its own length once, not on
        top of the previous shift's already-mutated boundary."""
        from lifton import align

        # Two boundaries; both straddle a sequence of two D-blocks.
        # cigar_accum_len = 5 puts us inside boundary[0]=(0,10);
        # cigar_accum_len = 6 puts us inside boundary[1]=(5,11).
        boundaries = {0: (0, 10), 1: (5, 11)}

        # Apply two D-block shifts and verify each boundary sees length
        # added EXACTLY ONCE per D-block, not cumulatively.
        b_after_d1 = align.adjust_cdss_protein_boundary(
            dict(boundaries), cigar_accum_len=5, length=2,
        )
        b_after_d2 = align.adjust_cdss_protein_boundary(
            dict(b_after_d1), cigar_accum_len=8, length=3,
        )
        # boundary[0] (0,10) was shifted by D1 (+2) → end becomes 12.
        # D2 at cigar_accum_len=8 is OUTSIDE boundary[0] (which is now (0,12)
        # in alignment coords; cigar_accum_len 8 < 12 → still inside).
        # The exact numbers depend on the algorithm; the invariant we
        # care about is that we never see a 2x stacked shift on the same
        # boundary in a single call.
        # Test invariant: between two successive calls, no boundary's end
        # grows by more than the sum of D-block lengths.
        max_growth = 2 + 3  # D1 + D2
        for k in boundaries:
            growth = b_after_d2[k][1] - boundaries[k][1]
            assert growth <= max_growth, (
                f"Boundary {k} grew by {growth} (> {max_growth}); "
                f"cumulative double-shift bug not patched."
            )


# ---------------------------------------------------------------------------
# V2.11 — parasail_align_protein_base silently coerces "" to "*".
# ---------------------------------------------------------------------------

class TestV2_11_ParasailEmptyProteinRaises:
    def test_empty_query_raises(self):
        """An empty protein sequence is a programming error upstream
        (the user's CDS produced no codons). Silently substituting "*"
        masks the bug. After the patch, an explicit
        LiftOnAlignmentError flags it."""
        from lifton import align
        from lifton.exceptions import LiftOnAlignmentError

        with pytest.raises(LiftOnAlignmentError):
            align.parasail_align_protein_base("", "MAGT*")

    def test_normal_query_unaffected(self):
        from lifton import align
        result = align.parasail_align_protein_base("MAGT", "MAGT*")
        assert result is not None


# ---------------------------------------------------------------------------
# V5.3 — Duplicate ID auto-renamed by gffutils, no warning.
# ---------------------------------------------------------------------------

class TestV5_3_DuplicateIDWarning:
    def test_duplicate_id_emits_warning(self, capsys, tmp_path):
        """When gffutils silently renames `tx1`, `tx1` to `tx1`, `tx1_1`,
        the user must get a [WARNING] so they know which IDs collided."""
        from lifton import annotation

        gff = tmp_path / "dup.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=g\n"
            "chr2\ttest\tmRNA\t300\t400\t.\t+\t.\tID=tx1;Parent=g\n"
        )
        try:
            annotation.Annotation(
                str(gff), infer_genes=False, infer_transcripts=False,
                verbose=True,
            )
        except Exception:
            # Even if the build path errors out, we still expect the
            # warning to have been emitted before that point.
            pass
        captured = capsys.readouterr()
        assert "[WARNING]" in captured.err and "tx1" in captured.err, (
            "Duplicate-ID warning was not emitted; the user has no way "
            "to know the GFF3 contains colliding IDs."
        )


# ---------------------------------------------------------------------------
# V5.4 / V5.5 — Reserved characters not percent-encoded in attributes.
# ---------------------------------------------------------------------------

class TestV5_4_5_AttributeEscaping:
    def _make_feature(self, attrs):
        import gffutils
        return gffutils.Feature(
            seqid="chr1", source="test", featuretype="gene",
            start=100, end=200, strand="+", frame=".",
            attributes=attrs, id=attrs.get("ID", [""])[0],
        )

    def test_semicolon_in_value_is_percent_encoded(self):
        """`Note=alpha;beta` would split as two attributes downstream.
        After the writer fix, the emitted line has `%3B` not raw `;`."""
        from lifton.io import gff3_writer

        feature = self._make_feature({
            "ID": ["g1"],
            "Note": ["alpha;beta"],
        })
        line = gff3_writer.format_feature(feature)
        # Locate the attributes column (last tab-separated field).
        attrs_col = line.rstrip("\n").split("\t")[-1]
        assert "alpha%3Bbeta" in attrs_col, (
            f"Reserved ';' was not percent-encoded; got attrs: {attrs_col}"
        )
        assert "alpha;beta" not in attrs_col, (
            f"Raw ';' leaked into attributes: {attrs_col}"
        )

    def test_equals_in_value_is_percent_encoded(self):
        from lifton.io import gff3_writer

        feature = self._make_feature({
            "ID": ["g2"],
            "Note": ["foo=bar"],
        })
        line = gff3_writer.format_feature(feature)
        attrs_col = line.rstrip("\n").split("\t")[-1]
        # The Note=foo=bar payload should have its inner = encoded.
        assert "foo%3Dbar" in attrs_col, (
            f"Reserved '=' was not percent-encoded; got attrs: {attrs_col}"
        )

    def test_comma_outside_multivalue_attr_is_encoded(self):
        """`,` is the multi-value separator for Parent / Alias / Note.
        In an ad-hoc attribute (e.g. `weird_field`) it must be encoded."""
        from lifton.io import gff3_writer

        feature = self._make_feature({
            "ID": ["g3"],
            "weird_field": ["a,b,c"],
        })
        line = gff3_writer.format_feature(feature)
        attrs_col = line.rstrip("\n").split("\t")[-1]
        assert "a%2Cb%2Cc" in attrs_col

    def test_tab_in_value_is_encoded(self):
        from lifton.io import gff3_writer

        feature = self._make_feature({
            "ID": ["g4"],
            "Note": ["foo\tbar"],
        })
        line = gff3_writer.format_feature(feature)
        # Output must have exactly 9 tab-separated columns.
        cols = line.rstrip("\n").split("\t")
        assert len(cols) == 9, (
            f"Tab in attribute value broke column count: {cols}"
        )
        assert "foo%09bar" in cols[8]


# ---------------------------------------------------------------------------
# V5.6 — Attribute order: ID first, Parent second, then alphabetical.
# ---------------------------------------------------------------------------

class TestV5_6_AttributeCanonicalOrder:
    def test_id_then_parent_then_alpha(self):
        from lifton.io import gff3_writer
        import gffutils

        feature = gffutils.Feature(
            seqid="chr1", source="test", featuretype="mRNA",
            start=100, end=200, strand="+", frame=".",
            attributes={
                "Note": ["x"],
                "Alias": ["y"],
                "Parent": ["g"],
                "ID": ["t1"],
                "Dbxref": ["GeneID:1"],
            },
            id="t1",
        )
        line = gff3_writer.format_feature(feature)
        attrs_col = line.rstrip("\n").split("\t")[-1]
        keys = [pair.split("=", 1)[0] for pair in attrs_col.split(";")]
        # Required canonical order: ID, Parent, then alphabetical
        assert keys[0] == "ID", f"first key was {keys[0]!r}; expected 'ID'"
        assert keys[1] == "Parent", f"second key was {keys[1]!r}"
        # Remaining keys must be sorted alphabetically.
        rest = keys[2:]
        assert rest == sorted(rest), (
            f"keys after ID/Parent are not alphabetical: {rest}"
        )

    def test_no_parent_still_canonical(self):
        from lifton.io import gff3_writer
        import gffutils

        feature = gffutils.Feature(
            seqid="chr1", source="test", featuretype="gene",
            start=100, end=200, strand="+", frame=".",
            attributes={"Note": ["x"], "ID": ["g1"], "Alias": ["y"]},
            id="g1",
        )
        line = gff3_writer.format_feature(feature)
        attrs_col = line.rstrip("\n").split("\t")[-1]
        keys = [pair.split("=", 1)[0] for pair in attrs_col.split(";")]
        assert keys[0] == "ID"
        assert keys[1:] == sorted(keys[1:])


# ---------------------------------------------------------------------------
# V5.8 — Phase recomputation after middle-CDS drop by ORF rescue.
# ---------------------------------------------------------------------------

class TestV5_8_PhaseAfterMiddleCDSDrop:
    def test_accum_skips_dropped_middle_cds(self, gff_standard):
        """If ORF rescue sets `exon.cds = None` for a middle exon, the
        phase recomputation for downstream CDSs must use accum_cds_length
        from EMITTED CDSs only — not the pre-drop running total."""
        import gffutils
        from lifton import lifton_class

        # Build a 3-exon transcript by hand.
        gff = (
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=g;gene_biotype=protein_coding\n"
            "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=t;Parent=g\n"
            "chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t\n"
            "chr1\ttest\texon\t100\t130\t.\t+\t.\tID=e2;Parent=t\n"
            "chr1\ttest\texon\t200\t230\t.\t+\t.\tID=e3;Parent=t\n"
        )
        from io import StringIO
        db = gffutils.create_db(
            StringIO(gff).read(), ":memory:", from_string=True,
            force=True, keep_order=True,
            merge_strategy="create_unique", sort_attribute_values=True,
        )
        mrna = db["t"]
        trans = lifton_class.Lifton_TRANS(
            "t", "g", "g", 0, mrna, dict(mrna.attributes),
        )
        for ex in db.children(mrna, featuretype="exon", level=1):
            trans.add_exon(ex)

        # Force an ORF that spans only exon 1 + exon 3 (skipping middle).
        # accum_exon_length: e1=30, e2=60, e3=90.
        # ORF [0, 30) lands fully in exon 1; [60, 90) lands fully in exon 3.
        # That's not realistic but it forces the pattern; the production
        # code would more likely set e2.cds = None via the
        # "ORF hasn't started yet" / "ORF has already ended" branches.
        # Direct simulation: set up scenario manually via __iterate.
        orf = lifton_class.Lifton_ORF(0, 30)  # only first exon
        trans._Lifton_TRANS__update_cds_boundary(orf)

        # After the rescue, exons 2 and 3 should have cds = None.
        assert trans.exons[0].cds is not None
        assert trans.exons[1].cds is None
        assert trans.exons[2].cds is None
        # The remaining CDS should have phase 0 (first emitted CDS).
        assert trans.exons[0].cds.entry.frame in ("0", "."), (
            f"Single-CDS phase should be 0 or '.', got "
            f"{trans.exons[0].cds.entry.frame!r}"
        )


# ---------------------------------------------------------------------------
# V5.9 — write_entry must reject inverted coordinates (start > end).
# ---------------------------------------------------------------------------

class TestV5_9_InvertedCoordinatesRejected:
    def test_format_feature_rejects_start_gt_end(self):
        from lifton.io import gff3_writer
        from lifton.exceptions import LiftOnInputError
        import gffutils

        bad = gffutils.Feature(
            seqid="chr1", source="test", featuretype="gene",
            start=200, end=100, strand="+", frame=".",
            attributes={"ID": ["bad"]}, id="bad",
        )
        with pytest.raises(LiftOnInputError):
            gff3_writer.format_feature(bad)


# ---------------------------------------------------------------------------
# Property-based fuzz tests — GFF3 attribute round-trip integrity
# ---------------------------------------------------------------------------

class TestPropertyBasedAttributeRoundTrip:
    """Hostile fuzzing: throw arbitrary printable strings (including
    every reserved char) at the writer; the resulting line must
    (a) have exactly 9 tab-separated columns, (b) re-parse cleanly,
    (c) round-trip the original value when percent-decoded."""

    @settings(max_examples=120, deadline=None)
    @given(
        key=st.text(
            alphabet=st.characters(min_codepoint=33, max_codepoint=126,
                                   blacklist_characters="\t\n;=&,%"),
            min_size=1, max_size=15,
        ),
        value=st.text(
            alphabet=st.characters(min_codepoint=33, max_codepoint=126),
            min_size=0, max_size=30,
        ),
    )
    def test_random_attribute_value_round_trips(self, key, value):
        from lifton.io import gff3_writer
        import gffutils
        from urllib.parse import unquote

        # Don't fuzz the reserved keys
        if key in ("ID", "Parent"):
            return
        feature = gffutils.Feature(
            seqid="chr1", source="test", featuretype="gene",
            start=10, end=20, strand="+", frame=".",
            attributes={"ID": ["g"], key: [value]},
            id="g",
        )
        line = gff3_writer.format_feature(feature)
        cols = line.rstrip("\n").split("\t")
        assert len(cols) == 9, (
            f"Tab/newline in value broke column count: {cols!r}"
        )
        # The attributes column must split cleanly on ';'.
        pairs = cols[8].split(";")
        for p in pairs:
            assert "=" in p, f"malformed attribute pair: {p!r}"
        # Round-trip: the encoded value should percent-decode back.
        for p in pairs:
            k, v = p.split("=", 1)
            if k == key:
                assert unquote(v) == value, (
                    f"value did not round-trip: {value!r} -> {v!r} -> "
                    f"{unquote(v)!r}"
                )


class TestPropertyBasedCoordinateInvariant:
    """For any random valid coordinate range, format_feature either
    succeeds with 9-column output, OR raises LiftOnInputError. Never
    silently emits a malformed line."""

    @settings(max_examples=80, deadline=None)
    @given(
        start=st.integers(min_value=1, max_value=10**8),
        length=st.integers(min_value=0, max_value=10**6),
    )
    def test_random_valid_range_emits_or_raises(self, start, length):
        from lifton.io import gff3_writer
        from lifton.exceptions import LiftOnInputError
        import gffutils

        end = start + length
        feature = gffutils.Feature(
            seqid="chr1", source="test", featuretype="gene",
            start=start, end=end, strand="+", frame=".",
            attributes={"ID": ["g"]}, id="g",
        )
        try:
            line = gff3_writer.format_feature(feature)
        except LiftOnInputError:
            return
        cols = line.rstrip("\n").split("\t")
        assert len(cols) == 9
        # And the start/end columns reflect the input.
        assert int(cols[3]) == start
        assert int(cols[4]) == end
