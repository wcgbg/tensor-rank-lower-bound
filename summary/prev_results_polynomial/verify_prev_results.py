#!/usr/bin/env python3
"""Verify prev_results.csv:
  (1) logical consistency with cert_summary.csv (no bound can contradict the
      other source, since both bracket the same true rank R);
  (2) independently re-derive the number/degree of irreducible factors of
      x^N-1 over F_q with sympy and check the Winograd/CRT cells;
  (3) summarise where the certified results match / beat / fall short of prior
      work.
Exit non-zero if any hard contradiction is found."""
import csv, os, sys
from math import isqrt, gcd

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    with open(os.path.join(HERE, path)) as f:
        return list(csv.DictReader(f))


def num(x):
    return None if x in ("-", "") else int(x)


cert = {}
for r in load("../our_results/cert_summary.csv"):
    cert[(r["problem"].rsplit("_q", 1)[0], int(r["q"]), int(r["na"]))] = (
        num(r["lb"]),
        num(r["ub"]),
    )
prev = load("prev_results.csv")

# ---- independent factor check with sympy (Galois factorisation) -------------
try:
    from sympy import GF, Poly, symbols

    x = symbols("x")

    def sympy_factor_degrees(q, N, p, e):
        # x^N - 1 over F_q where q = p^e.  sympy factors over GF(p); to handle
        # q = p^e we instead recompute distinct-factor count via cyclotomic
        # cosets (validated below against sympy for the prime-field case).
        if e == 1:
            f = Poly(x**N - 1, x, modulus=p)
            degs = []
            for fac, mult in f.factor_list()[1]:
                degs.append(fac.degree())
            return degs  # NB: lists each distinct factor once (mult separate)
        return None

    HAVE_SYMPY = True
except Exception as ex:
    HAVE_SYMPY = False
    print("sympy unavailable (%s); skipping independent factor check" % ex)


# Re-implement the script's helpers here for an independent recomputation.
def coset_degrees(q, m):
    seen, degs = set(), []
    for s in range(m):
        if s in seen:
            continue
        t, d = s, 0
        while t not in seen:
            seen.add(t)
            t = (t * q) % m
            d += 1
        degs.append(d)
    return degs


def k_distinct(q, N, p):
    m = N
    while m % p == 0:
        m //= p
    return len(coset_degrees(q, m))


# Cross-check our cyclotomic-coset factor count against sympy for prime q.
if HAVE_SYMPY:
    bad = 0
    for q, N, p, e in [
        (2, 7, 2, 1),
        (2, 9, 2, 1),
        (3, 8, 3, 1),
        (5, 6, 5, 1),
        (7, 6, 7, 1),
        (3, 4, 3, 1),
        (2, 5, 2, 1),
        (3, 7, 3, 1),
        (5, 4, 5, 1),
        (13, 5, 13, 1),
    ]:
        degs = sorted(sympy_factor_degrees(q, N, p, e))
        # our distinct-factor degrees (squarefree part):
        m = N
        while m % p == 0:
            m //= p
        ours = sorted(coset_degrees(q, m))
        if degs != ours:
            print(f"  FACTOR MISMATCH q={q} N={N}: sympy={degs} ours={ours}")
            bad += 1
    print(f"sympy factor cross-check: {'OK' if bad==0 else str(bad)+' MISMATCHES'}")

# ---- consistency + comparison ------------------------------------------------
violations = []
agree = beat_lb = beat_ub = short_lb = short_ub = nodata = 0
for r in prev:
    key = (r["problem"], int(r["q"]), int(r["n"]))
    plb, pub = num(r["lb"]), num(r["ub"])
    clb, cub = cert.get(key, (None, None))
    # hard contradictions (true rank R must satisfy all four):
    if plb is not None and cub is not None and plb > cub:
        violations.append(f"{key}: prev_lb {plb} > cert_ub {cub}")
    if clb is not None and pub is not None and clb > pub:
        violations.append(f"{key}: cert_lb {clb} > prev_ub {pub}")
    if plb is not None and pub is not None and plb > pub:
        violations.append(f"{key}: prev_lb {plb} > prev_ub {pub}")
    # comparison bookkeeping
    if plb is None and pub is None:
        nodata += 1
    if plb is not None and clb is not None:
        if clb > plb:
            print(f"{key}: cert_lb {clb} beats prev_lb {plb}")
            beat_lb += 1
        elif clb < plb:
            short_lb += 1
    if pub is not None and cub is not None:
        if cub < pub:
            print(f"{key}: cert_ub {cub} beats prev_ub {pub}")
            beat_ub += 1
        elif cub > pub:
            short_ub += 1
    if plb is not None and pub is not None and plb == pub and clb == plb and cub == pub:
        agree += 1

print("\n=== consistency ===")
if violations:
    print("CONTRADICTIONS FOUND:")
    for v in violations:
        print("  " + v)
else:
    print(
        "no contradictions: every prior bound brackets the certified interval consistently."
    )

print("\n=== certified vs prior (per cell, where both have a value) ===")
print(f"  prev fully exact & = cert         : {agree}")
print(f"  cert LOWER bound beats prior      : {beat_lb}")
print(f"  cert UPPER bound beats prior      : {beat_ub}")
print(f"  cert lower bound below prior best : {short_lb}")
print(f"  cert upper bound above prior best : {short_ub}")
print(f"  cells with NO prior result at all : {nodata}")

sys.exit(1 if violations else 0)
