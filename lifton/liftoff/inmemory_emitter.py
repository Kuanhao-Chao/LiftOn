# ---------------------------------------------------------------------------
# Phase 8 — in-memory GFF3 / GTF emitter for the vendored Liftoff.
#
# This is the pure-Python counterpart of write_new_gff.write_new_gff:
# same logic, same byte-for-byte output, but writes to a buffer in RAM
# instead of a file on disk. The Phase 8 streaming pipeline calls this
# to avoid the SQLite re-ingest round-trip after Liftoff finishes.
# ---------------------------------------------------------------------------

from __future__ import annotations

import io
import sys
from typing import Optional

from lifton.liftoff import liftoff_utils, write_new_gff
from lifton import __version__


def _write_header_to(buf, out_type: str, sys_argv: Optional[list] = None) -> None:
    """Mirrors :func:`write_new_gff.write_header` exactly."""
    if out_type == "gff3":
        buf.write("##gff-version 3" + "\n")
    buf.write("# LiftOn v" + __version__ + "\n")
    if sys_argv is None:
        sys_argv = sys.argv
    buf.write("# " + " ".join(sys_argv) + "\n")


def lifted_features_to_gff3_bytes(lifted_features, args, feature_db,
                                  *, sys_argv: Optional[list] = None) -> bytes:
    """Serialise the in-memory ``lifted_feature_list`` to a UTF-8 byte
    string in the same byte layout that :func:`write_new_gff.write_new_gff`
    would have written to disk.

    Used by the Phase 8 ``--inmemory-liftoff`` path so the parent
    pipeline can build a :class:`gffbase.FeatureDB` directly from the
    bytes without ever touching the disk.

    Parameters
    ----------
    lifted_features : dict
        ``{copy_id: [child Feature, ...]}`` — comes straight from
        :func:`lifton.liftoff.liftoff_main.run_all_liftoff_steps_inmemory`.
    args : argparse.Namespace
        Liftoff arguments. ``args.a`` and ``args.s`` thresholds are
        consulted by ``finalize_parent_features``.
    feature_db : gffutils.FeatureDB
        Liftoff's reference DB; ``feature_db.dialect['fmt']`` chooses
        ``gff3`` vs ``gtf`` output, exactly mirroring the disk path.
    sys_argv : list, optional
        Override for the ``# <command line>`` provenance comment.
        Defaults to ``sys.argv`` so the emitted bytes match what the
        on-disk file would look like for the same invocation.

    Returns
    -------
    bytes
        UTF-8-encoded GFF3 (or GTF) blob, identical byte-for-byte to
        the file ``write_new_gff`` would have written.
    """
    buf = io.StringIO()
    out_type = feature_db.dialect["fmt"]
    _write_header_to(buf, out_type, sys_argv=sys_argv)

    parents = liftoff_utils.get_parent_list(lifted_features)
    parents.sort(key=lambda x: x.id)
    final_parent_list = write_new_gff.finalize_parent_features(parents, args)
    final_parent_list.sort(key=lambda x: (x.seqid, x.start))

    for final_parent in final_parent_list:
        child_features = lifted_features[final_parent.attributes["copy_id"][0]]
        parent_child_dict = write_new_gff.build_parent_dict(child_features, final_parent)
        write_new_gff.write_feature(
            [final_parent], buf, child_features, parent_child_dict, out_type,
        )

    return buf.getvalue().encode("utf-8")
