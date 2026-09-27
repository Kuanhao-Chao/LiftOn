"""Count what a run drops, so the number reaches the user.

A lift can fail to emit a feature for many ordinary reasons: a reference
transcript the annotation does not actually declare, a sequence that cannot be
extracted, a database row that is missing where the code expected one. Each of
those sites already handled its own failure sensibly -- skip the feature, keep
going -- and most of them logged a line. What none of them did was *count*.

That distinction is not academic here. The Liftoff ``-copies`` resolution bug
dropped a transcript and all of its exons and CDS at one such site, ~4,400
times across the 17-genome benchmark corpus, in every release up to v1.0.11.
The warning was printed every single time -- 550 lines of it in one rice run --
and nobody saw it, because a warning per event in a stream of hundreds of
thousands of lines is indistinguishable from noise. A single number at the end
is not.

``locus_pipeline._record_childless_gene`` added exactly that for one class.
This module generalises it so every class of loss is counted the same way, in
one place, and reported together.

State is module-level rather than threaded through call signatures: several of
these sites sit deep in functions that have no ``args`` to hang a counter on
(``extract_sequence._stream_inner``), and adding a parameter to each of them
purely to carry a tally would be worse than a module global. ``reset()`` is
called once at the start of every run, which matters because the test suite
runs many pipelines in a single process.

Step 7 dispatches to a *thread* pool, so several of these sites are reached
concurrently and the tally is taken under a lock. A forked worker has its own
copy of this module, so the isoform scorer's workers keep a journal
(:func:`start_journal`) that returns with their results and is merged here
(:func:`merge`): the same scoring runs in-process at ``-t 1``, and the counts
must not depend on the thread count.

A feature is counted once per class, however many passes revisit it; a record
without an id is counted every time.
"""
from __future__ import annotations

import threading

#: Ids kept per class, so a user can start debugging without re-running.
EXAMPLE_CAP = 10

#: Every class of loss, with the sentence the summary prints. A ``record()``
#: for a kind that is not declared here raises, which is deliberate: a counter
#: nobody can interpret is barely better than no counter.
CLASSES: dict[str, str] = {
    "unresolvable_transcript":
        "reference transcript could not be resolved (the transcript and all "
        "its exons and CDS were dropped)",
    "unresolvable_gene":
        "reference gene id could not be resolved (the whole gene was dropped)",
    "reference_transcript_sequence":
        "reference transcript sequence could not be extracted (no DNA evidence "
        "for this transcript)",
    "reference_protein_sequence":
        "reference protein could not be extracted (no protein evidence, so "
        "miniprot cannot support or rescue this gene)",
    "miniprot_reference_missing":
        "reference gene or transcript row was missing while preparing a "
        "miniprot candidate (the candidate was scored without it)",
    "cross_locus_candidate":
        "cross-locus rescue candidate could not be read (it cannot carry its "
        "isoforms)",
    "hierarchy_depth_exceeded":
        "feature nested deeper than the pre-fetch depth limit (its own row is "
        "emitted, but every transcript, exon and CDS below it was dropped)",
    "miniprot_hit_unmapped":
        "miniprot hit could not be tied to any reference transcript (its "
        "Target= names a protein the id map does not contain)",
    "miniprot_gene_unresolved":
        "miniprot hit named a known reference transcript, but no reference "
        "gene could be found for it, so the hit could not be rescued",
    "reference_feature_length_missing":
        "reference gene's CDS span was not in the length index, so the "
        "length checks could not be applied and the candidate was skipped",
    "rescue_candidate_error":
        "miniprot-only rescue candidate raised while being built or scored "
        "(the error was recorded as a pipeline failure and the candidate "
        "abandoned)",
    "cds_spanning_exons":
        "CDS spanned more than one exon and could not be split at exon "
        "boundaries unambiguously, so the model carrying it was rejected "
        "rather than written across an intron or with coding bases counted "
        "twice",
}

_lock = threading.Lock()
_counts: dict[str, int] = {}
_examples: dict[str, list[str]] = {}
_seen: dict[str, set] = {}
#: Records kept by a forked worker for its parent to merge; None elsewhere.
_journal: list | None = None


def reset() -> None:
    """Start a fresh tally. Called once per run."""
    with _lock:
        _counts.clear()
        _examples.clear()
        _seen.clear()


def record(kind: str, feature_id=None) -> None:
    """Tally one dropped feature."""
    if kind not in CLASSES:
        raise KeyError(
            f"unknown drop class {kind!r}; declare it in drop_ledger.CLASSES "
            f"so the summary can explain it")
    with _lock:
        if _journal is not None:
            _journal.append((kind, feature_id))
        if feature_id is not None:
            seen = _seen.setdefault(kind, set())
            if str(feature_id) in seen:
                return
            seen.add(str(feature_id))
        _counts[kind] = _counts.get(kind, 0) + 1
        if feature_id is not None:
            kept = _examples.setdefault(kind, [])
            if len(kept) < EXAMPLE_CAP:
                kept.append(str(feature_id))


def start_journal() -> None:
    """Keep every record for a parent process to :func:`merge` (forked worker)."""
    global _journal
    with _lock:
        _journal = []


def take_journal() -> list:
    """The records kept since the last call, as ``[(kind, feature_id), ...]``."""
    global _journal
    with _lock:
        if _journal is None:
            return []
        taken, _journal = _journal, []
        return taken


def merge(entries) -> None:
    """Record what a forked worker recorded (duplicates count once)."""
    for kind, feature_id in entries or ():
        record(kind, feature_id)


def counts() -> dict[str, int]:
    """``{kind: n}`` for the classes that fired, in declaration order."""
    with _lock:
        return {kind: _counts[kind] for kind in CLASSES if _counts.get(kind)}


def examples(kind: str) -> list[str]:
    with _lock:
        return list(_examples.get(kind, ()))


def total() -> int:
    with _lock:
        return sum(_counts.values())


def report(manifest=None) -> None:
    """Print one summary block and record the totals in the run manifest.

    Silent when nothing was dropped -- a clean run should not grow a section
    telling the user that nothing happened.
    """
    tallied = counts()
    if manifest is not None:
        for kind in CLASSES:
            manifest.record_count(f"dropped_{kind}", _counts.get(kind, 0))
        manifest.record_count("dropped_features_total", total())
    if not tallied:
        return
    from lifton import logger
    lines = []
    for kind, count in tallied.items():
        shown = ", ".join(examples(kind)[:3])
        lines.append(f"{count} x {CLASSES[kind]}"
                     + (f"; e.g. {shown}" if shown else ""))
    logger.log_section(
        f"{total()} reference feature(s) were dropped", lines, kind="warning")
