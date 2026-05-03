# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""``GFFWriter`` — port of legacy ``gffutils.gffwriter.GFFWriter``."""

from __future__ import annotations

import os
import shutil
import tempfile
from typing import Iterable, Optional, Union

from .feature import Feature


class GFFWriter:
    """Write Feature records back to a GFF/GTF file."""

    def __init__(
        self,
        out: Union[str, "os.PathLike", "io.IOBase"],
        with_header: bool = True,
        in_place: bool = False,
    ):
        self.with_header = with_header
        self.in_place = in_place
        self._opened_path: Optional[str] = None
        self._target_path: Optional[str] = None

        if hasattr(out, "write"):
            self._fh = out
        elif in_place:
            # Atomic write via tempfile, swap on close.
            self._target_path = str(out)
            tmp = tempfile.NamedTemporaryFile(
                mode="w", delete=False, dir=os.path.dirname(self._target_path) or ".",
                suffix=".gffbase.tmp", encoding="utf-8",
            )
            self._fh = tmp.file
            self._opened_path = tmp.name
            tmp.close()
            self._fh = open(self._opened_path, "w", encoding="utf-8")
        else:
            self._opened_path = str(out)
            self._fh = open(self._opened_path, "w", encoding="utf-8")

        if self.with_header:
            self._fh.write("##gff-version 3\n")

    def write_rec(self, rec) -> None:
        if isinstance(rec, str):
            line = rec.rstrip("\n")
        else:
            line = str(rec)
        self._fh.write(line + "\n")

    def write_recs(self, recs: Iterable) -> None:
        for r in recs:
            self.write_rec(r)

    def write_gene_recs(self, db, gene_id) -> None:
        gene = db[gene_id] if isinstance(gene_id, str) else gene_id
        self.write_rec(gene)
        for child in db.children(gene, level=None, order_by="start"):
            self.write_rec(child)

    def write_mRNA_children(self, db, mrna_id) -> None:
        mrna = db[mrna_id] if isinstance(mrna_id, str) else mrna_id
        self.write_rec(mrna)
        for child in db.children(mrna, level=1, order_by="start"):
            self.write_rec(child)

    def write_exon_children(self, db, exon_id) -> None:
        exon = db[exon_id] if isinstance(exon_id, str) else exon_id
        self.write_rec(exon)
        for child in db.children(exon, level=1, order_by="start"):
            self.write_rec(child)

    def close(self) -> None:
        if self._fh is not None and not self._fh.closed:
            self._fh.flush()
            self._fh.close()
        if self.in_place and self._opened_path and self._target_path:
            shutil.move(self._opened_path, self._target_path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
