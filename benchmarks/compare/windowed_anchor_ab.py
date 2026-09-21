"""Is `_unique_anchors` cheaper when the reference index is restricted to the
k-mers the query actually contributes?

An anchor must be unique in BOTH sequences, so a k-mer the query does not
contribute can never become one and indexing it is pure cost. `candidate()`
below is the shipped implementation's shape; `current()` calls whatever the
tree currently has. Every case asserts the two return the SAME anchors before
either is timed -- a faster function that moves a single anchor would move a
window, and the windowed aligner is exact only because the windows are.

Run it from a checkout that predates the change to reproduce the "before"
column; against the current tree both arms are the same code and the ratio
collapses to 1, which is itself a check that the right thing shipped.
"""
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
from lifton import windowed_align


def current(query, ref, k):
    return windowed_align._unique_anchors(query, ref, k)


def candidate(query, ref, k):
    if len(query) < k or len(ref) < k:
        return []
    q_first, q_dup = {}, set()
    for i in range(len(query) - k + 1):
        km = query[i:i + k]
        if km in q_first:
            q_dup.add(km)
        else:
            q_first[km] = i
    for km in q_dup:
        del q_first[km]
    if not q_first:
        return []
    r_first, r_dup = {}, set()
    for i in range(len(ref) - k + 1):
        km = ref[i:i + k]
        if km in q_first:
            if km in r_first:
                r_dup.add(km)
            else:
                r_first[km] = i
    anchors = [(q_first[km], r) for km, r in r_first.items() if km not in r_dup]
    anchors.sort()
    return anchors


def mutate(seq, rate, alphabet, rng):
    out = []
    for ch in seq:
        r = rng.random()
        if r < rate * 0.7:
            out.append(rng.choice(alphabet))
        elif r < rate * 0.85:
            continue
        elif r < rate:
            out.append(ch); out.append(rng.choice(alphabet))
        else:
            out.append(ch)
    return "".join(out)


def case(name, alphabet, length, k, divergence, rng):
    ref = "".join(rng.choice(alphabet) for _ in range(length))
    query = mutate(ref, divergence, alphabet, rng)
    a = current(query, ref, k)
    b = candidate(query, ref, k)
    assert a == b, f"{name}: anchors differ ({len(a)} vs {len(b)})"
    reps = max(3, int(2_000_000 / length))
    t_cur = t_new = 0.0
    for _ in range(reps):                      # alternate, to share any drift
        s = time.perf_counter(); current(query, ref, k); t_cur += time.perf_counter() - s
        s = time.perf_counter(); candidate(query, ref, k); t_new += time.perf_counter() - s
    print(f"{name:28s} len={length:7d} k={k:2d} div={divergence:.2f} "
          f"anchors={len(a):6d}  current={t_cur/reps*1000:8.3f} ms  "
          f"candidate={t_new/reps*1000:8.3f} ms  {t_cur/t_new:.2f}x")


rng = random.Random(0)
AA = "ACDEFGHIKLMNPQRSTVWY"
NT = "ACGT"
print("protein (k=8):")
for length in (3000, 8000, 20000):
    for div in (0.02, 0.15):
        case("protein", AA, length, 8, div, rng)
print("\nDNA (k=15):")
for length in (10000, 50000, 120000):
    for div in (0.02, 0.15):
        case("dna", NT, length, 15, div, rng)
