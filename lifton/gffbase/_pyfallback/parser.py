# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Pure-Python streaming parser. Mirrors the Rust crate's behavior so it can
serve as a correctness oracle and as a fallback when the native extension is
unavailable.

Phase 16: same NCBI GFF3 validation rules as the Rust parser, with the same
`strict=True` raise / `strict=False` warning-collect semantics.
"""

from __future__ import annotations

import gzip
import io
from typing import Iterator, List, Optional

from ..dialect import merge_dialects
from ..feature import ParsedFeature
from .attributes import parse_attributes


def _gff_format_error_class():
    """Resolve the canonical ``GFFFormatError`` lazily.

    Importing it at module load creates a circular import (``gffbase``'s
    ``__init__`` re-binds the Python class to the Rust class when the
    extension is built — but that re-binding hasn't happened yet at the
    time this submodule loads). Resolving it on first use lets callers
    catch the *same* class regardless of whether the Rust extension or
    the Python fallback raised.
    """
    try:
        from .._native import GFFFormatError as _C   # type: ignore
        return _C
    except Exception:
        from ..exceptions import GFFFormatError as _C
        return _C


def _make_error(message: str, line_no: int, kind: str):
    cls = _gff_format_error_class()
    try:
        return cls(message, line_no=line_no, kind=kind)
    except TypeError:
        # The PyO3-derived class doesn't accept kwargs; build with the
        # positional message and attach the metadata as attributes
        # afterwards (matches what `gff_error_to_py` does in lib.rs).
        err = cls(message)
        try:
            err.line_no = line_no
            err.kind = kind
            err.message = message
        except (AttributeError, TypeError):
            pass
        return err


# Backwards-compat: tests importing `GFFFormatError` from this module.
class _LazyGFFFormatErrorProxy:
    def __instancecheck__(self, instance):
        return isinstance(instance, _gff_format_error_class())
    def __call__(self, *a, **kw):
        return _gff_format_error_class()(*a, **kw)


GFFFormatError = _gff_format_error_class()  # resolve eagerly enough for raise sites below


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _iter_lines(stream) -> Iterator[str]:
    for line in stream:
        if line.endswith("\n"):
            line = line[:-1]
        if line.endswith("\r"):
            line = line[:-1]
        yield line


def _parse_coord(s: str) -> Optional[int]:
    if s == "." or s == "":
        return None
    return int(s)


def _validate(
    *,
    line_no: int,
    seqid: str,
    featuretype: str,
    start: Optional[int],
    end: Optional[int],
    score: str,
    strand: str,
    frame: str,
    n_pairs: int,
    blob: str,
) -> Optional[GFFFormatError]:
    """Mirror of `validate.rs::validate_fields` + `validate_attributes_pairs`."""
    if not seqid:
        return _make_error(
            f"line {line_no}: seqid (col 1) is empty", line_no, "EmptySeqid",
        )
    if not featuretype:
        return _make_error(
            f"line {line_no}: featuretype (col 3) is empty",
            line_no, "EmptyFeaturetype",
        )
    if any(ch.isspace() for ch in featuretype):
        return _make_error(
            f"line {line_no}: featuretype contains whitespace: {featuretype!r}",
            line_no, "InvalidFeaturetype",
        )
    if start is not None and start < 1:
        return _make_error(
            f"line {line_no}: start coordinate must be >= 1 (got {start})",
            line_no, "InvalidCoordinate",
        )
    if start is not None and end is not None and end < start:
        return _make_error(
            f"line {line_no}: end < start ({end} < {start})",
            line_no, "InvalidCoordinate",
        )
    if strand not in ("+", "-", "?", "."):
        return _make_error(
            f"line {line_no}: strand must be one of '+', '-', '?', '.'; got {strand!r}",
            line_no, "InvalidStrand",
        )
    if frame not in (".", "0", "1", "2"):
        return _make_error(
            f"line {line_no}: phase must be 0, 1, 2, or '.'; got {frame!r}",
            line_no, "InvalidPhase",
        )
    if featuretype == "CDS" and frame == ".":
        return _make_error(
            f"line {line_no}: CDS row missing required phase (must be 0, 1, or 2)",
            line_no, "InvalidPhase",
        )
    if score not in ("", "."):
        try:
            float(score)
        except ValueError:
            return _make_error(
                f"line {line_no}: score must be a float or '.'; got {score!r}",
                line_no, "InvalidScore",
            )
    trimmed = blob.strip()
    if trimmed and trimmed != ".":
        # GFF3 attribute pairs use `key=value`; GTF uses `key "value"`.
        # The attribute parser is permissive enough to return a stub
        # `(token, "")` pair from raw garbage, so we require the blob
        # itself to show structure: at least one `=` (GFF3) OR a quoted
        # value (GTF) must appear, AND the parser must have produced
        # at least one pair.
        has_eq = "=" in trimmed
        has_quote = '"' in trimmed
        if (not has_eq and not has_quote) or n_pairs == 0:
            return _make_error(
                f"line {line_no}: attribute string did not parse into any "
                f"key=value pair: {trimmed[:60]!r}",
                line_no, "InvalidAttribute",
            )
    return None


def _coord_or_error(s: str, line_no: int, which: str) -> Optional[int]:
    """Returns the int, or raises GFFFormatError with structured info."""
    if s == "." or s == "":
        return None
    try:
        return int(s)
    except ValueError:
        raise _make_error(
            f"line {line_no}: {which} coordinate is not an integer: {s!r}",
            line_no, "InvalidCoordinate",
        )


def _parse_line_into_feature(line: str, line_no: int) -> ParsedFeature:
    fields = line.split("\t")
    if len(fields) < 9:
        raise _make_error(
            f"line {line_no}: expected at least 9 tab-separated fields, "
            f"found {len(fields)}",
            line_no, "TooFewFields",
        )
    seqid, source, featuretype, start_s, end_s, score, strand, frame = fields[:8]
    blob = fields[8]
    extra = fields[9:]
    pairs, _obs = parse_attributes(blob)
    start = _coord_or_error(start_s, line_no, "start")
    end = _coord_or_error(end_s, line_no, "end")
    err = _validate(
        line_no=line_no,
        seqid=seqid, featuretype=featuretype,
        start=start, end=end,
        score=score, strand=strand, frame=frame,
        n_pairs=len(pairs), blob=blob,
    )
    if err is not None:
        raise err
    return ParsedFeature(
        seqid=seqid,
        source=source,
        featuretype=featuretype,
        start=start,
        end=end,
        score=score,
        strand=strand,
        frame=frame,
        attributes_blob=blob.encode("utf-8"),
        attributes_pairs=pairs,
        extra=extra,
    )


def _stream_features(
    stream,
    checklines: int,
    force_dialect_check: bool,
    force_gff: bool,
    strict: bool,
    warnings: List[dict],
    directives: List[str],
):
    """Two-pass iteration: collect the first `checklines` features and their
    dialect observations, then continue streaming. Behaves identically when the
    file is shorter than `directives`. ``directives`` is mutated in place so
    callers can read it even if the file contains zero feature rows."""
    samples: List[dict] = []
    buffered: List[ParsedFeature] = []
    fasta_reached = False
    line_no = 0

    def _maybe_handle(err) -> bool:
        """Return True when caller should skip the line, False when caller
        should propagate (i.e., raise)."""
        if strict:
            return False
        warnings.append({
            "line_no": getattr(err, "line_no", 0),
            "kind":    getattr(err, "kind", ""),
            "message": getattr(err, "message", str(err)),
        })
        return True

    for line in _iter_lines(stream):
        line_no += 1
        if not line:
            continue
        if line.startswith("##"):
            if line.startswith("##FASTA"):
                fasta_reached = True
                break
            directives.append(line)
            continue
        if line.startswith("#"):
            continue
        try:
            feat = _parse_line_into_feature(line, line_no)
        except _gff_format_error_class() as e:
            if _maybe_handle(e):
                continue
            raise
        _, obs = parse_attributes(feat.attributes_blob.decode("utf-8", errors="replace"))
        samples.append(obs)
        buffered.append(feat)
        if not force_dialect_check and len(buffered) >= checklines:
            break

    dialect = merge_dialects(samples)
    if force_gff:
        dialect["fmt"] = "gff3"
        dialect["keyval separator"] = "="

    for f in buffered:
        yield f, directives, dialect

    if fasta_reached:
        return

    for line in _iter_lines(stream):
        line_no += 1
        if not line:
            continue
        if line.startswith("##"):
            if line.startswith("##FASTA"):
                return
            directives.append(line)
            continue
        if line.startswith("#"):
            continue
        try:
            feat = _parse_line_into_feature(line, line_no)
        except _gff_format_error_class() as e:
            if _maybe_handle(e):
                continue
            raise
        yield feat, directives, dialect


class _FallbackIterator:
    """Mirrors the Rust iterator's surface: __iter__/__next__, dialect(),
    directives(). Backwards-compat callers expect both forms.

    Phase 16: also exposes ``warnings`` (a list of dicts populated when
    ``strict=False``).
    """

    def __init__(self, stream, checklines: int, force_dialect_check: bool,
                 force_gff: bool, strict: bool = True):
        self._warnings: List[dict] = []
        self._directives: List[str] = []
        self._gen = _stream_features(
            stream, checklines, force_dialect_check, force_gff,
            strict=strict, warnings=self._warnings,
            directives=self._directives,
        )
        self._dialect = None
        self._exhausted = False

    def __iter__(self):
        return self

    def __next__(self) -> ParsedFeature:
        feat, directives, dialect = next(self._gen)
        self._dialect = dialect
        return feat

    def _drain_for_metadata(self) -> None:
        """Force the generator to its first yield (or to exhaustion) so
        that ``directives`` and ``dialect`` are populated even on
        directives-only / empty inputs. Idempotent."""
        if self._exhausted:
            return
        try:
            # The generator runs the dialect-peek phase before its first
            # yield, populating directives + dialect along the way. We
            # drive it just enough to reach that point.
            first = next(self._gen)
            # Restore: stash the first record so __next__ still sees it.
            self._gen = self._chain([first], self._gen)
        except StopIteration:
            self._exhausted = True

    @staticmethod
    def _chain(prefix, suffix):
        for item in prefix:
            yield item
        for item in suffix:
            yield item

    def dialect(self) -> dict:
        if self._dialect is None:
            self._drain_for_metadata()
        return self._dialect or {}

    def directives(self) -> List[str]:
        if not self._directives:
            self._drain_for_metadata()
        return list(self._directives)

    @property
    def warnings(self) -> List[dict]:
        return list(self._warnings)


def parse_file(
    path: str,
    checklines: int = 10,
    force_dialect_check: bool = False,
    force_gff: bool = False,
    strict: bool = True,
) -> _FallbackIterator:
    stream = _open(path)
    return _FallbackIterator(stream, checklines, force_dialect_check, force_gff, strict)


def parse_bytes(
    data: bytes,
    checklines: int = 10,
    force_dialect_check: bool = False,
    force_gff: bool = False,
    strict: bool = True,
) -> _FallbackIterator:
    stream = io.StringIO(data.decode("utf-8", errors="replace"))
    return _FallbackIterator(stream, checklines, force_dialect_check, force_gff, strict)


def detect_dialect(path: str, checklines: int = 10) -> dict:
    # Dialect detection is non-strict by design.
    it = parse_file(path, checklines=checklines, strict=False)
    drained: List[ParsedFeature] = []
    try:
        for _ in range(checklines):
            drained.append(next(it))
    except StopIteration:
        pass
    return it.dialect()
