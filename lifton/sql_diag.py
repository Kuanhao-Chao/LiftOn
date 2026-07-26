"""Attribute SQLite statements to the LiftOn call site that issued them.

``sqlite3.Cursor.execute`` is the third-heaviest entry in both Step-7 profiles
(24.0 s on drosophila, 15.1 s on dog->cat) after the parasail kernel and
``copy.deepcopy``. The profile says *how much* but not *which query* or *from
where*, and the per-locus materialiser issues several different ``children()``
shapes per feature. Guessing which one to batch is exactly the mistake the
Step-7 profile already corrected twice, so measure first.

Entirely opt-in and output-neutral: nothing is installed unless
``LIFTON_SQL_DIAG=<path>`` is set, and the trace callback only observes.

    LIFTON_SQL_DIAG=/tmp/sql.tsv lifton ... -t 1

The report groups by (normalised SQL, nearest ``lifton/`` frame), so a hot
statement points straight at the code that has to change.
"""
from __future__ import annotations

import collections
import os
import re
import sys
import threading
import traceback


_LOCK = threading.Lock()
_COUNTS: "collections.Counter[tuple]" = collections.Counter()
_INSTALLED: set = set()

# Collapse literals so the same query shape aggregates.
_LITERAL = re.compile(r"'[^']*'|\b\d+\b")
_WS = re.compile(r"\s+")


def enabled() -> bool:
    return bool(os.environ.get("LIFTON_SQL_DIAG"))


def _normalise(statement: str) -> str:
    return _WS.sub(" ", _LITERAL.sub("?", statement)).strip()[:180]


def _call_site() -> str:
    """The innermost ``lifton/`` frame that is not this module."""
    for frame in reversed(traceback.extract_stack()):
        name = frame.filename
        if "/lifton/" in name and not name.endswith("sql_diag.py"):
            short = name.split("/lifton/", 1)[1]
            return f"{short}:{frame.lineno} {frame.name}"
    return "<outside lifton>"


def attach(db, label: str) -> None:
    """Start recording statements issued through ``db``'s connection."""
    if not enabled():
        return
    conn = getattr(db, "conn", None)
    if conn is None or not hasattr(conn, "set_trace_callback"):
        return
    if id(conn) in _INSTALLED:
        return
    _INSTALLED.add(id(conn))

    def _record(statement: str) -> None:
        key = (label, _normalise(statement or ""), _call_site())
        with _LOCK:
            _COUNTS[key] += 1

    try:
        conn.set_trace_callback(_record)
    except Exception:                                  # pragma: no cover
        _INSTALLED.discard(id(conn))


def report() -> None:
    """Write the tally to ``LIFTON_SQL_DIAG`` and summarise on stderr."""
    path = os.environ.get("LIFTON_SQL_DIAG")
    if not path:
        return
    with _LOCK:
        rows = _COUNTS.most_common()
    total = sum(count for _, count in rows)
    try:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("count\tshare\tdb\tcall_site\tstatement\n")
            for (label, statement, site), count in rows:
                handle.write(
                    f"{count}\t{count / (total or 1):.4f}\t{label}\t{site}\t"
                    f"{statement}\n"
                )
    except OSError as exc:                              # pragma: no cover
        sys.stderr.write(f"[LiftOn][sql] could not write {path}: {exc}\n")
        return
    sys.stderr.write(
        f"\n[LiftOn][sql] {total} statements over {len(rows)} "
        f"(statement, call-site) pairs -> {path}\n"
    )
    for (label, statement, site), count in rows[:12]:
        sys.stderr.write(
            f"[LiftOn][sql]   {count:9d} {count / (total or 1):6.1%}  "
            f"{label:9s} {site:46s} {statement[:70]}\n"
        )
    sys.stderr.flush()
