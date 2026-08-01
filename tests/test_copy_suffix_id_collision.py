"""A repeated id must not be renamed onto a feature that already exists.

gffutils' ``create_unique`` disambiguates a repeated ``cds-X`` by renaming it to
``cds-X_1``, ``cds-X_2``, … But ``_<n>`` is exactly the suffix Liftoff's
``-copies`` mode gives extra gene copies, so on a copies-bearing annotation that
generated key can already belong to a different feature. The insert then dies
with ``UNIQUE constraint failed: features.id`` and the whole run aborts.

v1.0.8 survived this because its fallback transform pre-renamed repeats into a
``_dup<n>`` namespace that cannot collide. That transform was later reduced to a
no-op to stop rewriting logical ids -- a good goal, but it left the fallback with
nothing to fall back to, and three plant genomes that v1.0.8 lifts stopped
building at all. Found by re-benchmarking the release: v1.0.8 built 34 of 34
benchmark inputs, v1.0.10 built 31.

The shape below is the real one, reduced: a CDS in two segments sharing
``ID=cds-A`` (spec-valid discontinuous CDS), plus a copy feature that already
owns ``cds-A_1`` -- the exact name gffutils would generate.

Order is load-bearing. gffutils renames on insert, so the copy feature has to be
inserted *before* the repeat it collides with; with the copy last, gffutils
renames the repeat to a name that is still free and nothing breaks. Liftoff
emits copies ahead of the primary locus, which is why the real files trip it.
"""

from __future__ import annotations

import gffutils
import pytest

from lifton import annotation


COLLIDING_GFF3 = """\
##gff-version 3
chr1\tLiftOn\tgene\t8000\t9000\t.\t+\t.\tID=gene-A_1;extra_copy_number=1
chr1\tLiftOn\tmRNA\t8000\t9000\t.\t+\t.\tID=rna-A_1;Parent=gene-A_1
chr1\tLiftOn\texon\t8000\t9000\t.\t+\t.\tID=exon-A_1-1;Parent=rna-A_1
chr1\tLiftOn\tCDS\t8000\t9000\t.\t+\t0\tID=cds-A_1;Parent=rna-A_1
chr1\tLiftOn\tgene\t1000\t3000\t.\t+\t.\tID=gene-A
chr1\tLiftOn\tmRNA\t1000\t3000\t.\t+\t.\tID=rna-A;Parent=gene-A
chr1\tLiftOn\texon\t1000\t1500\t.\t+\t.\tID=exon-A-1;Parent=rna-A
chr1\tLiftOn\texon\t2500\t3000\t.\t+\t.\tID=exon-A-2;Parent=rna-A
chr1\tLiftOn\tCDS\t1000\t1500\t.\t+\t0\tID=cds-A;Parent=rna-A
chr1\tLiftOn\tCDS\t2500\t3000\t.\t+\t0\tID=cds-A;Parent=rna-A
"""

# Same discontinuous CDS, but nothing occupies `cds-A_1`: gffutils can rename
# on its own and LiftOn must not interfere.
CLEAN_GFF3 = """\
##gff-version 3
chr1\tLiftOn\tgene\t1000\t3000\t.\t+\t.\tID=gene-B
chr1\tLiftOn\tmRNA\t1000\t3000\t.\t+\t.\tID=rna-B;Parent=gene-B
chr1\tLiftOn\texon\t1000\t1500\t.\t+\t.\tID=exon-B-1;Parent=rna-B
chr1\tLiftOn\texon\t2500\t3000\t.\t+\t.\tID=exon-B-2;Parent=rna-B
chr1\tLiftOn\tCDS\t1000\t1500\t.\t+\t0\tID=cds-B;Parent=rna-B
chr1\tLiftOn\tCDS\t2500\t3000\t.\t+\t0\tID=cds-B;Parent=rna-B
"""


def _write(tmp_path, text, name="in.gff3"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_plain_gffutils_cannot_ingest_the_colliding_shape(tmp_path):
    """Establish that the fixture really does defeat gffutils unaided.

    Without this the test below could pass for the wrong reason.
    """
    path = _write(tmp_path, COLLIDING_GFF3)
    with pytest.raises(Exception, match="UNIQUE constraint failed"):
        gffutils.create_db(
            str(path), str(tmp_path / "plain.db"),
            merge_strategy="create_unique", force=True, verbose=False,
        )


def test_lifton_ingests_the_colliding_shape(tmp_path):
    """The regression: this aborted the run on released v1.0.10."""
    path = _write(tmp_path, COLLIDING_GFF3)
    ann = annotation.Annotation(str(path), None, False)
    assert ann.db_connection is not None


def test_the_preexisting_copy_feature_survives(tmp_path):
    """`cds-A_1` is a real feature and must keep its own identity."""
    path = _write(tmp_path, COLLIDING_GFF3)
    ann = annotation.Annotation(str(path), None, False)
    ids = {feature.id for feature in ann.db_connection.all_features()}
    assert "cds-A_1" in ids
    # The repeats are renamed into a namespace that cannot collide.
    assert not any(i.startswith("cds-A_1_dup") for i in ids)


def test_discontinuous_cds_is_left_alone_when_nothing_collides(tmp_path):
    """Only contaminated ids get rewritten; ordinary shared ids do not.

    This is the property the no-op transform was reaching for, and the reason
    the fix renames selectively instead of restoring the old blanket rename.
    """
    path = _write(tmp_path, CLEAN_GFF3)
    ann = annotation.Annotation(str(path), None, False)
    ids = {feature.id for feature in ann.db_connection.all_features()}
    assert not any("_dup" in i for i in ids), (
        "a shared CDS id with no colliding copy feature must not be renamed"
    )


def test_contamination_is_detected_from_the_existing_scan(tmp_path):
    """The predicate reuses the preflight scan rather than re-reading the file."""
    path = _write(tmp_path, COLLIDING_GFF3)
    ann = annotation.Annotation(str(path), None, False)
    assert ann._contaminated_ids() == {"cds-A"}

    clean = _write(tmp_path, CLEAN_GFF3, name="clean.gff3")
    ann_clean = annotation.Annotation(str(clean), None, False)
    assert ann_clean._contaminated_ids() == set()


def test_a_bug_in_the_fallback_is_not_reported_as_a_bad_file(tmp_path):
    """A coding error must surface, not masquerade as a DB-build failure.

    The strategy handlers catch `Exception` so a malformed input falls through
    to the next strategy. That breadth also swallowed a NameError raised inside
    the transform itself while this fix was being written, which showed up as
    "DB build failed" and cost real debugging time.
    """
    path = _write(tmp_path, COLLIDING_GFF3)
    ann = annotation.Annotation.__new__(annotation.Annotation)

    with pytest.raises(NameError):
        ann._reraise_if_our_bug(NameError("undefined helper"))
    # Genuine database errors still fall through to the next strategy.
    assert ann._reraise_if_our_bug(ValueError("malformed record")) is None
