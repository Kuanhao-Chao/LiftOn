# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""``DataIterator`` — legacy-compatible factory wrapping the new parser.

Yields ``Feature`` objects (not raw ``ParsedFeature``) so downstream code
that consumes the iterator and prints features Just Works.
"""

from __future__ import annotations

from typing import Iterator, List, Optional

from . import parser as _parser
from .feature import Feature, ParsedFeature


class _DataIterator:
    """Public iterator. Exposes ``.dialect`` and ``.directives`` like legacy."""

    def __init__(
        self,
        data,
        checklines: int = 10,
        transform=None,
        force_dialect_check: bool = False,
        from_string: bool = False,
        **kwargs,
    ):
        if from_string:
            self._inner = _parser.parse_bytes(
                data.encode("utf-8") if isinstance(data, str) else data,
                checklines=checklines,
                force_dialect_check=force_dialect_check,
            )
        else:
            self._inner = _parser.parse_gff(
                data,
                checklines=checklines,
                force_dialect_check=force_dialect_check,
            )
        self._transform = transform

    def __iter__(self) -> Iterator[Feature]:
        return self

    def __next__(self) -> Feature:
        pf: ParsedFeature = next(self._inner)
        feat = Feature(
            seqid=pf.seqid,
            source=pf.source,
            featuretype=pf.featuretype,
            start=pf.start,
            end=pf.end,
            score=pf.score,
            strand=pf.strand,
            frame=pf.frame,
            attributes=pf.attributes_blob,
            extra=("\t".join(pf.extra)) if pf.extra else None,
            dialect=self._inner.dialect() or {"fmt": "gff3"},
        )
        if self._transform is not None:
            out = self._transform(feat)
            if out is False:
                return self.__next__()
            if out is not None:
                feat = out
        return feat

    @property
    def dialect(self) -> dict:
        return self._inner.dialect() or {"fmt": "gff3"}

    @property
    def directives(self) -> List[str]:
        return list(self._inner.directives())


def DataIterator(
    data,
    checklines: int = 10,
    transform=None,
    force_dialect_check: bool = False,
    from_string: bool = False,
    **kwargs,
) -> _DataIterator:
    """Legacy factory. Returns an iterator yielding ``Feature``."""
    return _DataIterator(
        data,
        checklines=checklines,
        transform=transform,
        force_dialect_check=force_dialect_check,
        from_string=from_string,
        **kwargs,
    )
