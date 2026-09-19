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


#: Ordinal -> replacement, for ``str.translate``. Same mapping the
#: per-character loop applied, expressed so CPython can walk the string in C.
#: Profiling Step 7 on whole genomes (notes/step7_profile_2026-09.md) counted
#: 5.8 M calls to this function on rice and 17.3 M on dog to cat -- 5 % and 3 %
#: of the phase -- and every one of them ran a Python-level loop over every
#: character of every attribute value, almost always to return it unchanged.
_ENCODE_TRANSLATION = {ord("%"): "%25"}
_ENCODE_TRANSLATION.update(
    {ord(ch): "%{:02X}".format(ord(ch)) for ch in RESERVED_CHARS}
)

#: 86 % of real attribute values contain nothing that needs encoding, so the
#: cheapest correct answer is to notice that in C and hand back the string
#: itself. Translating unconditionally is only 1.3x the old loop; checking
#: first is 3.2x.
_NEEDS_ENCODING = re.compile(
    "[" + re.escape("".join(sorted(RESERVED_CHARS | {"%"}))) + "]"
)


def encode_attribute_value(value: str) -> str:
    """Percent-encode reserved characters per NCBI GFF3 spec.

    The encoding is unambiguous: every literal `%` becomes `%25` so a
    decoder (`urllib.parse.unquote`) can recover the original input
    byte-exactly. This means a user-supplied value containing the
    literal text "%00" round-trips to "%00", not to a NUL byte.
    """
    if value is None:
        return ""
    s = str(value)
    if not s or _NEEDS_ENCODING.search(s) is None:
        return s
    return s.translate(_ENCODE_TRANSLATION)


# Backward-compatible private name retained for focused tests and in-tree
# callers written before the encoder became shared with vendored Liftoff.
_encode_reserved = encode_attribute_value


def _attr_value_str(values: Any) -> str:
    """Render an attribute's value(s) into the GFF3 wire format.

    A list is joined with `,` (per NCBI multi-value rule). Each
    individual value is percent-encoded so any embedded `;`, `=`, `,`,
    `&`, tab, or newline doesn't break the line structure.
    """
    if isinstance(values, (list, tuple)):
        encoded = [encode_attribute_value(str(v)) for v in values]
        return ",".join(encoded)
    return encode_attribute_value(str(values))


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


def format_directives(directives) -> str:
    """V5.7 / Phase 15a — render the directive block that prefixes the
    GFF3 output. Always begins with `##gff-version 3`; preserves any
    additional `##` / `#!` directives passed in (in input order),
    de-duplicating the gff-version line if it appears more than once.

    Returns a string ending with a single newline so the caller can
    `fw.write(format_directives(...))` directly before the first
    feature row.
    """
    seen: set[str] = set()
    out: list[str] = ["##gff-version 3"]
    seen.add("##gff-version 3")
    for raw in directives or ():
        if raw is None:
            continue
        line = str(raw).rstrip("\r\n")
        if not line:
            continue
        if not (line.startswith("##") or line.startswith("#!")):
            # Not a directive — silently skip; callers may pass mixed
            # input as a convenience.
            continue
        if line in seen:
            continue
        seen.add(line)
        out.append(line)
    return "\n".join(out) + "\n"


def target_output_directives(directives):
    """Return only directives that remain true after assembly lift-over.

    Coordinate ranges, species/build pragmas, annotation provenance, and a
    reference ``##FASTA`` marker describe the *input* assembly. Re-emitting
    them on target coordinates produces a misleading output header. Ontology
    declarations are assembly-independent and are therefore retained.
    """
    allowed_prefixes = (
        "##gff-version",
        "##feature-ontology",
        "##attribute-ontology",
        "##source-ontology",
    )
    return [
        str(raw).rstrip("\r\n")
        for raw in directives or ()
        if raw is not None
        and str(raw).rstrip("\r\n").startswith(allowed_prefixes)
    ]


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
