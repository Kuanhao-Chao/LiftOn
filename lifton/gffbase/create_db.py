# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""``create_db()`` — drop-in successor to ``gffutils.create_db``.

Wraps ``gffbase.ingest.from_file`` with the legacy signature so downstream
code that does ``gffutils.create_db(path, ":memory:")`` works unchanged.
Many legacy kwargs are accepted-but-no-op for now and will be wired up in
Phase 6 (e.g. ``merge_strategy``, ``id_spec``, ``transform``).
"""

from __future__ import annotations

import os
import tempfile
from typing import Optional

from . import ingest as _ingest
from .interface import FeatureDB


def create_db(
    data,
    dbfn,
    *,
    id_spec=None,
    force: bool = False,
    verbose: bool = False,
    checklines: int = 10,
    from_string: bool = False,
    force_gff: bool = False,
    force_dialect_check: bool = False,
    merge_strategy: str = "error",
    force_merge_fields=None,
    transform=None,
    gtf_transcript_key: str = "transcript_id",
    gtf_gene_key: str = "gene_id",
    gtf_subfeature: str = "exon",
    disable_infer_genes: bool = False,
    disable_infer_transcripts: bool = False,
    infer_gene_extent: bool = True,
    keep_order: bool = False,
    text_factory=str,
    pragmas: Optional[dict] = None,
    sort_attribute_values: bool = False,
    dialect: Optional[dict] = None,
    _keep_tempfiles: bool = False,
    **kwargs,
) -> FeatureDB:
    """Create a database from a GFF3/GTF source.

    For Phase 5 the actively-honored kwargs are: ``data``, ``dbfn``, ``force``,
    ``checklines``, ``from_string``, ``disable_infer_genes``,
    ``disable_infer_transcripts``, ``gtf_subfeature``. The rest are accepted
    for signature compatibility and will be wired up in Phase 6.
    """
    cleanup_path: Optional[str] = None
    if from_string:
        # Materialize to a temp file so the parser can mmap-style read it.
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".gff3", delete=False, encoding="utf-8"
        )
        tmp.write(data)
        tmp.close()
        path = tmp.name
        cleanup_path = path
    else:
        path = data

    try:
        con, stats = _ingest.from_file(
            path,
            dbfn=dbfn,
            force=force,
            disable_infer_genes=disable_infer_genes,
            disable_infer_transcripts=disable_infer_transcripts,
            gtf_subfeature=gtf_subfeature,
        )
    finally:
        if cleanup_path and not _keep_tempfiles:
            try:
                os.unlink(cleanup_path)
            except OSError:
                pass

    db = FeatureDB((con, stats), keep_order=keep_order,
                   sort_attribute_values=sort_attribute_values,
                   text_factory=text_factory, pragmas=pragmas)
    return db
