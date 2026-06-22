"""
windowed_align.py — anchor-windowed global alignment for very long sequences.

LiftOn aligns every transcript with `parasail.nw_trace_scan_sat`, a full
Needleman–Wunsch DP + traceback that is O(L²) in BOTH time and memory. On giant
genes this explodes: a single titin DNA alignment (~106 kb) is ~106k × 106k ≈
1.1e10 cells whose traceback matrix is tens of GB (measured: the mouse benchmark
peaked at ~48 GB and took 1264 s, almost all of it in a handful of giant genes —
see memory `lifton-giant-gene-align-blowup`). That is both a ~10× slowdown and an
OOM crash on any machine under ~48 GB.

This module computes a near-diagonal global alignment in O(L·band)-ish time and
bounded memory by:
  1. seeding with k-mers that are UNIQUE in both sequences (reliable anchors),
  2. chaining them into one co-linear, non-overlapping diagonal (LIS),
  3. aligning only the small windows BETWEEN anchors with full parasail and
     concatenating — anchors themselves contribute exact-match columns.

The result is returned as a small object that is drop-in compatible with the
parasail result LiftOn consumes: it exposes `.traceback.{query,comp,ref}` and a
`.cigar` whose `.decode` is the CIGAR bytes (so `res.cigar.decode.decode()`
works, matching parasail's `=`/`X`/`I`/`D` convention — `I` = gap in ref, `D` =
gap in query). The CIGAR is derived from the produced alignment strings, so it is
always self-consistent with `.traceback`.

For near-identical (same/cross-species) lifts this reproduces full-DP identity
essentially exactly (validated on real titin: Δ identity = +0.00000, 82× faster).
As of Iteration 3 "band everything" is the DEFAULT (gate 2500 aa / 8000 nt, cap
1500; see `lifton/align.py:configure_alignment`), so the whole O(L²) mid-tail is
windowed — exact wherever anchors exist or a divergent region is below the giant
boundary (`max_fulldp`). Below the gate (tiny genes / fixtures) the exact full-DP
path is unchanged. `--full-dp-align` restores the exact giant-only path.
"""

import os

import parasail

# Defaults. Unique k-mers: AA has a 20-letter alphabet so k=8 is unique within
# any realistic protein; DNA has 4 letters so it needs a longer seed (k=15 →
# 4^15 ≈ 1e9 >> any transcript length). MIN_K is the floor when a window is
# anchor-sparse and we shrink k to find finer anchors.
DEFAULT_K_AA = 8
DEFAULT_K_NT = 15
MIN_K_AA = 4
MIN_K_NT = 8
# A window this small (or smaller) is aligned directly with full parasail. Caps
# any single DP allocation at ~WINDOW_CAP² cells (~tens of MB), so peak memory
# stays bounded regardless of input length. Env-tunable so the "band everything"
# fast-align mode (lifton/align.py) can lower the band granularity in-process
# without touching the giant-only default.
WINDOW_CAP = int(os.environ.get("LIFTON_ALIGN_WINDOW_CAP", "6000"))
_MAX_DEPTH = 24


def set_window_cap(n):
    """Set the full-DP leaf size (band granularity) in-process. Used by
    `lifton.align.configure_alignment` to switch to a finer band without an
    env round-trip; functions read the module global, so a plain rebind here
    is picked up immediately."""
    global WINDOW_CAP
    WINDOW_CAP = int(n)


class _Cigar:
    """Mimics parasail's `result.cigar`: `.decode` is the CIGAR bytes."""
    __slots__ = ("decode",)

    def __init__(self, cigar_bytes):
        self.decode = cigar_bytes


class _Traceback:
    __slots__ = ("query", "comp", "ref")

    def __init__(self, query, comp, ref):
        self.query = query
        self.comp = comp
        self.ref = ref


class _ShimResult:
    """parasail-result-shaped: exposes `.traceback` and `.cigar`."""
    __slots__ = ("traceback", "cigar")

    def __init__(self, query_aln, comp, ref_aln):
        self.traceback = _Traceback(query_aln, comp, ref_aln)
        self.cigar = _Cigar(_cigar_from_aln(query_aln, ref_aln))


def _cigar_from_aln(query_aln, ref_aln):
    """Derive a parasail-convention CIGAR (bytes) from aligned strings.

    Symbols: `=` match, `X` mismatch, `I` gap in ref (query base present),
    `D` gap in query (ref base present) — matching parasail.nw_trace_scan_sat.
    """
    ops = []  # (symbol, run_length)
    for q, r in zip(query_aln, ref_aln):
        if q == "-":
            sym = "D"          # gap in query
        elif r == "-":
            sym = "I"          # gap in ref
        elif q == r:
            sym = "="
        else:
            sym = "X"
        if ops and ops[-1][0] == sym:
            ops[-1][1] += 1
        else:
            ops.append([sym, 1])
    return ("".join(f"{n}{s}" for s, n in ops)).encode()


def _unique_anchors(query, ref, k):
    """Return anchors (q_pos, r_pos) for k-mers occurring exactly once in BOTH
    sequences (so they are unambiguous), sorted by q_pos."""
    if len(query) < k or len(ref) < k:
        return []
    q_count, q_first = {}, {}
    for i in range(len(query) - k + 1):
        km = query[i:i + k]
        c = q_count.get(km, 0) + 1
        q_count[km] = c
        if c == 1:
            q_first[km] = i
    r_count, r_first = {}, {}
    for i in range(len(ref) - k + 1):
        km = ref[i:i + k]
        c = r_count.get(km, 0) + 1
        r_count[km] = c
        if c == 1:
            r_first[km] = i
    anchors = [(q_first[km], r_first[km])
               for km, c in q_count.items()
               if c == 1 and r_count.get(km) == 1]
    anchors.sort()
    return anchors


def _chain_colinear(anchors, k):
    """Longest co-linear, non-overlapping anchor chain.

    anchors are sorted by q (q strictly increasing). Take the longest strictly
    r-increasing subsequence (patience sorting), then drop anchors that overlap
    the previous kept one (spacing < k in either axis)."""
    if not anchors:
        return []
    import bisect
    tails_r, tails_i, prev = [], [], [-1] * len(anchors)
    for i, (_q, r) in enumerate(anchors):
        pos = bisect.bisect_left(tails_r, r)
        if pos == len(tails_r):
            tails_r.append(r)
            tails_i.append(i)
        else:
            tails_r[pos] = r
            tails_i[pos] = i
        prev[i] = tails_i[pos - 1] if pos > 0 else -1
    chain = []
    idx = tails_i[-1]
    while idx != -1:
        chain.append(anchors[idx])
        idx = prev[idx]
    chain.reverse()
    # enforce non-overlap in both axes
    kept = []
    last_q = last_r = -1
    for q, r in chain:
        if q >= last_q + k and r >= last_r + k:
            kept.append((q, r))
            last_q, last_r = q, r
    return kept


def _full(query, ref, gap_open, gap_extend, matrix):
    res = parasail.nw_trace_scan_sat(query, ref, gap_open, gap_extend, matrix)
    return res.traceback.query, res.traceback.comp, res.traceback.ref


def _blind_split(query, ref, gap_open, gap_extend, matrix, cap):
    """Last-resort bounded alignment when a large window has no usable anchors
    (anchor-sparse divergent stretch). Split query into <=cap chunks and ref
    proportionally, align each chunk pair with full parasail, concatenate.
    Approximate but memory-bounded; only reached for pathological inputs."""
    Q, C, R = [], [], []
    n = max(1, (len(query) + cap - 1) // cap)
    for j in range(n):
        qa, qb = j * len(query) // n, (j + 1) * len(query) // n
        ra, rb = j * len(ref) // n, (j + 1) * len(ref) // n
        qs, rs = query[qa:qb], ref[ra:rb]
        if not qs and not rs:
            continue
        if not qs:
            Q.append("-" * len(rs)); C.append(" " * len(rs)); R.append(rs)
        elif not rs:
            Q.append(qs); C.append(" " * len(qs)); R.append("-" * len(qs))
        else:
            q_a, c_a, r_a = _full(qs, rs, gap_open, gap_extend, matrix)
            Q.append(q_a); C.append(c_a); R.append(r_a)
    return "".join(Q), "".join(C), "".join(R)


def _align_region(query, ref, gap_open, gap_extend, matrix, k, min_k, depth,
                  max_fulldp):
    """Recursively align (query, ref). Returns (query_aln, comp, ref_aln).

    `max_fulldp` is the largest anchor-less region we are willing to align with
    exact full DP before falling back to the approximate (but memory-bounded)
    blind split. With the giant-only default it equals WINDOW_CAP, so the
    fallback never fires (anything reaching it is already > WINDOW_CAP) and the
    Iteration-2 path is unchanged; "band everything" raises it to the giant
    boundary so divergent non-giant regions stay exact."""
    if query == "" and ref == "":
        return "", "", ""
    if query == "":                      # pure insertion of ref → gap in query
        return "-" * len(ref), " " * len(ref), ref
    if ref == "":                        # pure deletion → gap in ref
        return query, " " * len(query), "-" * len(query)
    if max(len(query), len(ref)) <= WINDOW_CAP or depth >= _MAX_DEPTH:
        return _full(query, ref, gap_open, gap_extend, matrix)

    anchors = _chain_colinear(_unique_anchors(query, ref, k), k)
    if not anchors:
        if k > min_k:
            return _align_region(query, ref, gap_open, gap_extend, matrix,
                                 max(min_k, k // 2), min_k, depth + 1, max_fulldp)
        # No anchors even at the finest k. Prefer exact full DP for regions up to
        # the memory-safe bound; only true giants take the approximate split.
        if max(len(query), len(ref)) <= max_fulldp:
            return _full(query, ref, gap_open, gap_extend, matrix)
        return _blind_split(query, ref, gap_open, gap_extend, matrix, WINDOW_CAP)

    nxt_k = max(min_k, k // 2)
    Q, C, R = [], [], []
    pq = pr = 0
    for qi, ri in anchors:
        wq, wc, wr = _align_region(query[pq:qi], ref[pr:ri], gap_open,
                                   gap_extend, matrix, nxt_k, min_k, depth + 1,
                                   max_fulldp)
        Q.append(wq); C.append(wc); R.append(wr)
        # exact-match anchor columns (k-mer identical in both by construction)
        Q.append(query[qi:qi + k]); C.append("|" * k); R.append(ref[ri:ri + k])
        pq, pr = qi + k, ri + k
    wq, wc, wr = _align_region(query[pq:], ref[pr:], gap_open, gap_extend,
                               matrix, nxt_k, min_k, depth + 1, max_fulldp)
    Q.append(wq); C.append(wc); R.append(wr)
    return "".join(Q), "".join(C), "".join(R)


def windowed_traceback(query, ref, gap_open, gap_extend, matrix, *,
                       k=None, min_k=None, is_dna=False, max_fulldp=None):
    """Drop-in replacement for `parasail.nw_trace_scan_sat(query, ref, ...)` for
    long sequences. Returns a parasail-result-shaped object exposing
    `.traceback.{query,comp,ref}` and `.cigar.decode` (bytes).

    `max_fulldp` bounds the exact-DP fallback for anchor-less regions; when
    None it defaults to WINDOW_CAP (fallback off → giant-only behaviour)."""
    if k is None:
        k = DEFAULT_K_NT if is_dna else DEFAULT_K_AA
    if min_k is None:
        min_k = MIN_K_NT if is_dna else MIN_K_AA
    if max_fulldp is None:
        max_fulldp = WINDOW_CAP
    q_aln, comp, r_aln = _align_region(query, ref, gap_open, gap_extend,
                                       matrix, k, min_k, 0, max_fulldp)
    return _ShimResult(q_aln, comp, r_aln)
