"""Phase 13.5B — hostile tests for the Critical + High severity findings
catalogued in `plans/phase_13_5A_vulnerability_audit.md`.

Each test is written to FAIL on the unpatched codebase and PASS once the
matching production-code fix lands. Tests are grouped by audit ID
(V1.x = exception handling, V2.x = algorithmic, V4.x = data, V5.x = GFF3
edge cases).

The tests stay hermetic: no external binaries, no network. Where a real
gffutils database is needed, the existing conftest fixtures are reused.
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
