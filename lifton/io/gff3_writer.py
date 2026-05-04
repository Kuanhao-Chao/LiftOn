"""Canonical GFF3 line writer (Phase 13.5C).

One single source of truth for serialising a `gffutils.Feature` (or
anything quacking like it) into a GFF3 line. Replaces the inline
``str(self.entry)`` calls scattered across `Lifton_GENE / Lifton_TRANS
/ Lifton_EXON / Lifton_CDS / LiftOn_FEATURE.write_entry`.

Implements the NCBI § Attribute Specifications:
  * **Reserved-character percent-encoding** (V5.4 / V5.5)
    — `;`, `=`, `&`, `,`, `\\t`, `\\n`, `\\r`, plus `%` itself when not
      already part of a valid `%XX` escape.
  * **Canonical attribute order** (V5.6) — `ID`, `Parent`, then
    alphabetical by key.
  * **Coordinate invariant** (V5.9) — `start <= end`; raises
    `LiftOnInputError` on inversion.

The function returns the line WITHOUT a trailing newline so the
caller controls EOL bytes (downstream writers add `"\\n"`).
"""

from __future__ import annotations

import re
from typing import Any, Mapping

from lifton.exceptions import LiftOnInputError
from lifton.io.ncbi_gff3_spec import RESERVED_CHARS

# Already-encoded sequences look like %XX where X is hex.
_PCT_ESCAPE_RE = re.compile(r"%[0-9A-Fa-f]{2}")


def _encode_reserved(value: str) -> str:
    """Percent-encode reserved characters per NCBI GFF3 spec.

    The encoding is unambiguous: every literal `%` becomes `%25` so a
    decoder (`urllib.parse.unquote`) can recover the original input
    byte-exactly. This means a user-supplied value containing the
    literal text "%00" round-trips to "%00", not to a NUL byte.
    """
    if value is None:
        return ""
    s = str(value)
    if not s:
        return s
    out_chars: list[str] = []
    for ch in s:
        if ch == "%":
            out_chars.append("%25")
        elif ch in RESERVED_CHARS:
            out_chars.append("%{:02X}".format(ord(ch)))
        else:
            out_chars.append(ch)
    return "".join(out_chars)


def _attr_value_str(values: Any) -> str:
    """Render an attribute's value(s) into the GFF3 wire format.

    A list is joined with `,` (per NCBI multi-value rule). Each
    individual value is percent-encoded so any embedded `;`, `=`, `,`,
    `&`, tab, or newline doesn't break the line structure.
    """
    if isinstance(values, (list, tuple)):
        encoded = [_encode_reserved(str(v)) for v in values]
        return ",".join(encoded)
    return _encode_reserved(str(values))


def _canonical_attr_order(keys):
    """Sort attribute keys: ID first, Parent second, alphabetical
    after. Stable for keys appearing exactly once."""
    keys = list(keys)
    head: list[str] = []
    if "ID" in keys:
        head.append("ID")
        keys.remove("ID")
    if "Parent" in keys:
        head.append("Parent")
        keys.remove("Parent")
    return head + sorted(keys)


def format_attributes(attributes: Mapping[str, Any]) -> str:
    """Serialise an attribute dict into the column-9 string for a
    GFF3 row. Order and encoding follow NCBI canonical rules."""
    if not attributes:
        return ""
    parts: list[str] = []
    for key in _canonical_attr_order(attributes.keys()):
        val = attributes[key]
        # Skip empty values (e.g. None or empty list).
        if val is None:
            continue
        if isinstance(val, (list, tuple)) and not val:
            continue
        parts.append(f"{key}={_attr_value_str(val)}")
    return ";".join(parts)


def format_feature(feature) -> str:
    """Render a `gffutils.Feature`-like object as a single GFF3 line
    (without trailing newline). Validates the start/end invariant and
    canonicalises attribute order + escaping.
    """
    seqid = getattr(feature, "seqid", ".") or "."
    source = getattr(feature, "source", ".") or "."
    ftype = getattr(feature, "featuretype", ".") or "."
    start = getattr(feature, "start", None)
    end = getattr(feature, "end", None)
    score = getattr(feature, "score", ".") or "."
    strand = getattr(feature, "strand", ".") or "."
    frame = getattr(feature, "frame", ".") or "."

    # V5.9 — coordinate invariant
    if start is None or end is None:
        raise LiftOnInputError(
            f"format_feature: missing start/end on feature "
            f"{getattr(feature, 'id', '<unknown>')!r}"
        )
    try:
        s_int = int(start)
        e_int = int(end)
    except (TypeError, ValueError) as exc:
        raise LiftOnInputError(
            f"format_feature: non-integer coordinates on feature "
            f"{getattr(feature, 'id', '<unknown>')!r}: {exc}"
        )
    if s_int < 1:
        raise LiftOnInputError(
            f"format_feature: GFF3 coordinates must be >= 1; got "
            f"start={s_int} on feature "
            f"{getattr(feature, 'id', '<unknown>')!r}"
        )
    if s_int > e_int:
        raise LiftOnInputError(
            f"format_feature: start ({s_int}) > end ({e_int}) on "
            f"feature {getattr(feature, 'id', '<unknown>')!r} — "
            "violates NCBI § Column 4-5."
        )

    attrs_col = format_attributes(getattr(feature, "attributes", {}) or {})

    return "\t".join([
        str(seqid), str(source), str(ftype),
        str(s_int), str(e_int),
        str(score), str(strand), str(frame),
        attrs_col,
    ])
