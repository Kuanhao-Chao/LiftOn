# ---------------------------------------------------------------------------
# Phase 10 — common types for the native binding facade.
# ---------------------------------------------------------------------------

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass(frozen=True)
class MinimapHit:
    """Single minimap2 alignment hit, mirrored from mappy.Alignment.

    Carries only the columns that LiftOn (or a SAM round-trip) actually
    consumes — the source object exposes more, but those are the ones
    the differential parity tests assert against the subprocess path.
    """

    query_name: str
    ctg: str            # target contig name
    r_st: int           # 0-based reference start (mappy convention)
    r_en: int           # exclusive reference end
    q_st: int           # 0-based query start
    q_en: int           # exclusive query end
    strand: int         # +1 or -1
    mapq: int
    NM: int             # edit distance
    cigar_str: str      # CIGAR string
    is_primary: bool


@dataclass(frozen=True)
class GFF3Hit:
    """Single miniprot GFF3 line decoded into structured fields. The
    public shape matches what `pyminiprot` (when it exists) is expected
    to yield, AND what the existing miniprot subprocess output produces
    when parsed line-by-line. The facade can therefore swap from the
    subprocess path to a real PyO3 binding with no caller changes."""

    seqid: str
    source: str
    featuretype: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: str          # raw attribute string (e.g. "ID=MP1;Target=tx1 1 66")

    @classmethod
    def from_gff_line(cls, line: str) -> Optional["GFF3Hit"]:
        """Parse one tab-separated GFF3 row. Returns None for blank
        lines or comment/directive lines."""
        if not line or line.startswith("#"):
            return None
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            return None
        try:
            start = int(cols[3])
            end = int(cols[4])
        except ValueError:
            return None
        return cls(
            seqid=cols[0],
            source=cols[1],
            featuretype=cols[2],
            start=start,
            end=end,
            score=cols[5],
            strand=cols[6],
            phase=cols[7],
            attributes=cols[8],
        )


@dataclass
class GFF3Bundle:
    """Group of GFF3Hit records returned by MiniprotIndex.align(...).
    Carries the raw bytes blob too so downstream callers that prefer
    the existing gffbase ingest path can keep using it without
    reconstructing the original lines."""

    hits: List[GFF3Hit] = field(default_factory=list)
    raw_bytes: bytes = b""

    def __iter__(self):
        return iter(self.hits)

    def __len__(self):
        return len(self.hits)
