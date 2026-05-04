# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Feature objects.

This module hosts two distinct types:

* ``ParsedFeature`` — slotted dataclass emitted by the parser.
* ``Feature`` — full backward-compatible public class, mirroring legacy
  ``gffutils.Feature``. It supports lazy attribute parsing so that DB rows
  loaded but never inspected don't pay the cost of decoding col-9.
"""

from __future__ import annotations

import json
from collections.abc import MutableMapping
from dataclasses import dataclass, field
from typing import Iterator, List, Mapping, Optional, Tuple, Union


@dataclass(slots=True)
class ParsedFeature:
    seqid: str
    source: str
    featuretype: str
    start: Optional[int]
    end: Optional[int]
    score: str
    strand: str
    frame: str
    # Raw col-9 bytes preserved for byte-faithful round-trip in Phase 4.
    attributes_blob: bytes
    # Long-form (key, value, multivalue_index) triples. `idx` preserves
    # the order of multi-valued attributes (e.g., Parent=a,b,c becomes
    # three rows with idx 0,1,2).
    attributes_pairs: List[Tuple[str, str, int]] = field(default_factory=list)
    # Tab-separated columns past column 9, if any.
    extra: List[str] = field(default_factory=list)

    @property
    def chrom(self) -> str:
        return self.seqid

    @property
    def stop(self) -> Optional[int]:
        return self.end

    def attributes_dict(self) -> dict:
        """Materialize attributes as `{key: [values...]}`. Preserves first-seen
        key order and multi-value ordering. Defers to `attributes_pairs` so the
        Rust and Python parsers remain trivially comparable."""
        out: dict = {}
        for k, v, _idx in self.attributes_pairs:
            out.setdefault(k, []).append(v)
        return out

    @classmethod
    def from_tuple(cls, tup) -> "ParsedFeature":
        """Build from the 11-tuple shape that the Rust extension yields."""
        (
            seqid,
            source,
            featuretype,
            start,
            end,
            score,
            strand,
            frame,
            blob,
            pairs,
            extra,
        ) = tup
        return cls(
            seqid=seqid,
            source=source,
            featuretype=featuretype,
            start=start,
            end=end,
            score=score,
            strand=strand,
            frame=frame,
            attributes_blob=bytes(blob) if not isinstance(blob, bytes) else blob,
            attributes_pairs=[(k, v, int(i)) for (k, v, i) in pairs],
            extra=list(extra),
        )


# ---------------------------------------------------------------------------
# Public ``Feature`` class — backward-compatible with legacy ``gffutils.Feature``.
# ---------------------------------------------------------------------------


# Field order for integer-indexed __getitem__ / __setitem__ (legacy semantics).
_FIELD_ORDER = (
    "seqid", "source", "featuretype", "start", "end",
    "score", "strand", "frame", "attributes",
)


def _coord_to_int(v) -> Optional[int]:
    """Normalize a coordinate. Legacy accepts int, '.', '', or None."""
    if v is None or v == "" or v == ".":
        return None
    if isinstance(v, int):
        return v
    return int(v)


class _LazyAttributes(MutableMapping):
    """Dict-like attribute store. Values are always lists.

    If constructed with a raw col-9 ``blob``, parsing is deferred until first
    access. This honors the Phase 2 §3.3 invariant: attributes never decoded
    unless someone reads them.
    """

    __slots__ = ("_d", "_blob", "_dialect_fmt", "_parsed")

    def __init__(
        self,
        initial: Union[Mapping, List[Tuple[str, str, int]], None] = None,
        blob: Optional[bytes] = None,
        dialect_fmt: str = "gff3",
    ):
        self._d: dict = {}
        self._blob = blob
        self._dialect_fmt = dialect_fmt
        self._parsed = False
        if initial is not None:
            self._ingest(initial)
            self._parsed = True
            self._blob = None

    # ----- internal -----

    def _ingest(self, initial):
        if isinstance(initial, _LazyAttributes):
            initial._materialize()
            for k, v in initial._d.items():
                self._d[k] = list(v)
            return
        if isinstance(initial, Mapping):
            for k, v in initial.items():
                self._d[k] = v if isinstance(v, list) else [v]
            return
        if isinstance(initial, list):  # list of (key, value, idx)
            for triple in initial:
                if len(triple) == 3:
                    k, v, _idx = triple
                else:
                    k, v = triple
                self._d.setdefault(k, []).append(v)
            return
        raise TypeError(f"unsupported attributes init type: {type(initial)!r}")

    def _materialize(self):
        if self._parsed:
            return
        if self._blob is None or self._blob == b"":
            self._parsed = True
            return
        # Defer to the pure-Python attribute parser to avoid a Rust hop on a
        # single line. This is the warm path for users who *do* read attrs.
        from ._pyfallback.attributes import parse_attributes
        text = self._blob.decode("utf-8", errors="replace")
        pairs, _obs = parse_attributes(text)
        for k, v, _idx in pairs:
            self._d.setdefault(k, []).append(v)
        self._parsed = True
        self._blob = None

    # ----- MutableMapping -----

    def __getitem__(self, key: str):
        self._materialize()
        return self._d[key]

    def __setitem__(self, key: str, value):
        self._materialize()
        self._d[key] = value if isinstance(value, list) else [value]

    def __delitem__(self, key: str):
        self._materialize()
        del self._d[key]

    def __iter__(self) -> Iterator[str]:
        self._materialize()
        return iter(self._d)

    def __len__(self) -> int:
        self._materialize()
        return len(self._d)

    def __contains__(self, key) -> bool:
        self._materialize()
        return key in self._d

    def __repr__(self) -> str:
        self._materialize()
        return f"Attributes({self._d!r})"

    # ----- helpers -----

    def items(self):
        self._materialize()
        return self._d.items()

    def keys(self):
        self._materialize()
        return self._d.keys()

    def values(self):
        self._materialize()
        return self._d.values()


class Feature:
    """Backward-compatible public Feature object.

    Mirrors the legacy ``gffutils.Feature`` constructor and observable
    behavior: 1-based inclusive coordinates, list-wrapped multi-value
    attributes, dialect-faithful ``__str__`` round-trip.
    """

    __slots__ = (
        "seqid", "source", "featuretype",
        "start", "end", "score", "strand", "frame",
        "attributes", "extra",
        "bin", "id", "dialect", "file_order",
        "keep_order", "sort_attribute_values",
        "_attributes_blob",
        "children",  # populated by FeatureDB.merge to expose component features
    )

    def __init__(
        self,
        seqid: str = ".",
        source: str = ".",
        featuretype: str = ".",
        start=".",
        end=".",
        score: str = ".",
        strand: str = ".",
        frame: str = ".",
        attributes=None,
        extra=None,
        bin: Optional[int] = None,
        id: Optional[str] = None,
        dialect: Optional[dict] = None,
        file_order: Optional[int] = None,
        keep_order: bool = False,
        sort_attribute_values: bool = False,
    ):
        self.seqid = seqid
        self.source = source
        self.featuretype = featuretype
        self.start = _coord_to_int(start)
        self.end = _coord_to_int(end)
        self.score = score if score is not None else "."
        self.strand = strand if strand is not None else "."
        self.frame = frame if frame is not None else "."
        self.bin = bin
        self.id = id
        self.dialect = dialect or {}
        self.file_order = file_order
        self.keep_order = keep_order
        self.sort_attribute_values = sort_attribute_values
        self._attributes_blob = None
        self.children = None

        fmt = (self.dialect or {}).get("fmt", "gff3")
        if isinstance(attributes, _LazyAttributes):
            self.attributes = attributes
        elif isinstance(attributes, (bytes, bytearray)):
            self._attributes_blob = bytes(attributes)
            self.attributes = _LazyAttributes(blob=self._attributes_blob, dialect_fmt=fmt)
        elif attributes is None:
            self.attributes = _LazyAttributes(initial={}, dialect_fmt=fmt)
        else:
            self.attributes = _LazyAttributes(initial=attributes, dialect_fmt=fmt)

        if extra is None:
            self.extra = []
        elif isinstance(extra, (bytes, bytearray)):
            text = extra.decode("utf-8", errors="replace")
            self.extra = text.split("\t") if text else []
        elif isinstance(extra, str):
            self.extra = extra.split("\t") if extra else []
        else:
            self.extra = list(extra)

    # ----- aliases -----

    @property
    def chrom(self) -> str:
        return self.seqid

    @chrom.setter
    def chrom(self, v: str):
        self.seqid = v

    @property
    def stop(self) -> Optional[int]:
        return self.end

    @stop.setter
    def stop(self, v):
        self.end = _coord_to_int(v)

    # ----- dunders -----

    def __len__(self) -> int:
        if self.start is None or self.end is None:
            return 0
        return self.end - self.start + 1

    def __repr__(self) -> str:
        return (
            f"<Feature {self.featuretype} ({self.seqid}:{self.start}-{self.end}"
            f"[{self.strand}]) at {hex(id(self))}>"
        )

    def __str__(self) -> str:
        return self._format_line()

    def __unicode__(self) -> str:
        return self.__str__()

    def __hash__(self) -> int:
        return hash(self._format_line())

    def __eq__(self, other) -> bool:
        if not isinstance(other, Feature):
            return NotImplemented
        return self._format_line() == other._format_line()

    def __ne__(self, other) -> bool:
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        return not eq

    def __getitem__(self, key):
        if isinstance(key, int):
            return getattr(self, _FIELD_ORDER[key])
        return self.attributes[key]

    def __setitem__(self, key, value):
        if isinstance(key, int):
            setattr(self, _FIELD_ORDER[key], value)
        else:
            self.attributes[key] = value

    # ----- formatting -----

    def _format_attributes(self) -> str:
        # If we have the original col-9 bytes and the user hasn't materialized
        # / mutated the attributes mapping, re-emit them verbatim. This is the
        # byte-faithful round-trip path.
        if (
            self._attributes_blob is not None
            and isinstance(self.attributes, _LazyAttributes)
            and not self.attributes._parsed
        ):
            return self._attributes_blob.decode("utf-8", errors="replace")

        fmt = (self.dialect or {}).get("fmt", "gff3")
        sep = "; " if (self.dialect or {}).get("field separator") == "; " else ";"
        kv_sep = (self.dialect or {}).get("keyval separator") or ("=" if fmt == "gff3" else " ")
        multival = (self.dialect or {}).get("multival separator", ",")
        items = list(self.attributes.items())
        if self.sort_attribute_values:
            items = [(k, sorted(v)) for k, v in items]

        parts = []
        for k, vs in items:
            vs_list = vs if isinstance(vs, list) else [vs]
            if fmt == "gff3":
                joined = multival.join(str(v) for v in vs_list)
                parts.append(f"{k}{kv_sep}{joined}")
            else:
                # GTF: typically `key "value"; key "value";`
                quoted = (self.dialect or {}).get("quoted GFF2 values", True)
                for v in vs_list:
                    if quoted:
                        parts.append(f'{k}{kv_sep}"{v}"')
                    else:
                        parts.append(f"{k}{kv_sep}{v}")
        s = sep.join(parts)
        if (self.dialect or {}).get("trailing semicolon") and not s.endswith(";"):
            s = s + ";"
        if (self.dialect or {}).get("leading semicolon"):
            s = ";" + s
        return s

    def _format_line(self) -> str:
        start_s = "." if self.start is None else str(self.start)
        end_s = "." if self.end is None else str(self.end)
        cols = [
            self.seqid, self.source, self.featuretype,
            start_s, end_s, self.score, self.strand, self.frame,
            self._format_attributes(),
        ]
        cols.extend(self.extra)
        return "\t".join(cols)

    # ----- legacy methods -----

    def astuple(self, encoding=None):
        """Legacy 12-tuple shape used by the SQLite export path:
        ``(id, seqid, source, featuretype, start, end, score, strand, frame,
            attributes_json, extra_json, bin)``.
        """
        attrs_dict = {k: list(v) for k, v in self.attributes.items()}
        return (
            self.id,
            self.seqid,
            self.source,
            self.featuretype,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            json.dumps(attrs_dict, separators=(",", ":")),
            json.dumps(self.extra, separators=(",", ":")) if self.extra else "[]",
            self.bin if self.bin is not None else self.calc_bin(),
        )

    def calc_bin(self, _bin=None) -> Optional[int]:
        if _bin is not None:
            self.bin = _bin
            return _bin
        if self.start is None or self.end is None:
            return None
        from ._bins import bin_from_coords
        self.bin = bin_from_coords(self.start, self.end)
        return self.bin

    def sequence(self, fasta, use_strand: bool = True) -> str:
        """Extract sequence from a FASTA path or a pyfaidx-style mapping."""
        if isinstance(fasta, str):  # pragma: no cover - pyfaidx is optional
            try:
                import pyfaidx  # type: ignore
            except ImportError as e:
                raise ImportError(
                    "Feature.sequence(path=...) requires the optional `pyfaidx` package"
                ) from e
            fa = pyfaidx.Fasta(fasta)
        else:
            fa = fasta
        seq = str(fa[self.seqid][self.start - 1 : self.end])
        if use_strand and self.strand == "-":
            seq = _revcomp(seq)
        return seq


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Construction from a DuckDB row.
# ---------------------------------------------------------------------------

# Row column order produced by FeatureDB._yield_features.
_DB_ROW_FIELDS = (
    "id", "seqid", "source", "featuretype", "start", "end",
    "score", "strand", "frame", "attributes_blob", "extra_blob", "file_order",
)


def feature_from_row(row, dialect: Optional[dict] = None) -> Feature:
    """Build a ``Feature`` from a DuckDB row tuple. Lazy in attributes."""
    (
        fid, seqid, source, featuretype, start, end,
        score, strand, frame, blob, extra_blob, file_order,
    ) = row
    return Feature(
        seqid=seqid,
        source=source,
        featuretype=featuretype,
        start=start,
        end=end,
        score=score if score is not None else ".",
        strand=strand if strand is not None else ".",
        frame=frame if frame is not None else ".",
        attributes=blob if blob is not None else None,
        extra=extra_blob if extra_blob else None,
        id=fid,
        dialect=dialect or {"fmt": "gff3"},
        file_order=file_order,
    )
