"""LiftOn exception hierarchy (Phase 13.5B).

These types replace the bare `except:` and broad `except Exception:`
swallows that the Phase 13.5A audit flagged as Critical/High. Using a
named hierarchy lets callers catch *LiftOn* failures without also
swallowing KeyboardInterrupt, SystemExit, or unrelated programming
errors.
"""

from __future__ import annotations


class LiftOnError(Exception):
    """Root of the LiftOn exception hierarchy."""


class LiftOnInputError(LiftOnError):
    """Raised when the user-supplied annotation, FASTA, or auxiliary
    file is invalid or unprocessable. Always points at fixable user
    input — never at an internal logic bug."""


class LiftOnAlignmentError(LiftOnError):
    """Raised when sequence alignment cannot be carried out (e.g.
    parasail kernel failure, BioPython translation aborted, parsing
    of a CIGAR string fails)."""
