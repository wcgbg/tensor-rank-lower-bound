#!/usr/bin/env python3
"""Generate prev_results.csv: previously-published best-known lower/upper bounds
for the bilinear complexity (tensor rank) of the three multiplication problems
studied in this repo, for every (problem, p, m, q, n) cell that appears in
certs/cert_summary.csv.

    problem  : full | cyclic | extfield
    full     : product of two degree-(n-1) polynomials over F_q (2n-1 outputs)
               -> rank denoted  M_q(n)
    cyclic   : multiplication in F_q[x]/(x^n - 1)               (n outputs)
    extfield : multiplication in F_{q^n} over F_q               (n outputs)
               -> rank denoted  mu_q(n)

Design (as requested): *each reference in the literature is one Python function*.
A reference function inspects a cell and, if the theorem/lemma/table it encodes
says something about that cell, returns a Claim(kind, value, ref).  kind is
'lb' or 'ub'.  combine() then keeps, per cell, the largest lower bound and the
smallest upper bound, breaking ties toward the more specific reference.

Nothing here is invented: a cell is left blank ('-') when no prior publication is
known to the author of this script to bound it.  See REFERENCES at the bottom.
"""

import csv
import os
from math import isqrt, gcd
from collections import namedtuple

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "prev_results.csv")

# Grid to emit: all (problem, q=p^m, n) with these ranges.
PROBLEMS = ("cyclic", "extfield", "full")
Q_MAX = 16
N_MAX = 11

Claim = namedtuple("Claim", ["kind", "value", "ref", "prio"])
# prio: higher = more specific / preferred when two refs give the same value.

# --------------------------------------------------------------------------- #
#  Number-theoretic helpers                                                    #
# --------------------------------------------------------------------------- #


def prime_powers(hi):
    """List of (q, p, m) with q = p^m a prime power, lo <= q <= hi."""
    out = []
    for q in range(2, hi + 1):
        p = next(d for d in range(2, q + 1) if q % d == 0)  # smallest factor=prime
        t, m = q, 0
        while t % p == 0:
            t //= p
            m += 1
        if t == 1:  # q is a pure p^m
            out.append((q, p, m))
    return out


def is_square(q):
    r = isqrt(q)
    return r * r == q


def epsilon(q):
    """Shokrollahi's eps(q): if q is a square, 2*sqrt(q); else the greatest
    integer <= 2*sqrt(q) that is prime to q.  (Ballet survey, Thm 'thm_shokr'.)"""
    if is_square(q):
        return 2 * isqrt(q)
    t = isqrt(4 * q)  # floor(2*sqrt(q))
    while gcd(t, q) != 1:
        t -= 1
    return t


def cyclotomic_coset_degrees(q, m):
    """Degrees of the irreducible factors of x^m - 1 over F_q (gcd(m,q)=1):
    the sizes of the q-cyclotomic cosets modulo m.  Sum of degrees = m."""
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


def num_distinct_irreducible_factors(q, N, p):
    """k = number of *distinct* monic irreducible factors of x^N - 1 over F_q,
    where p = char(F_q).  Write N = p^a * mm with gcd(mm,p)=1; then
    x^N-1 = (x^mm-1)^(p^a) and k = #q-cyclotomic cosets mod mm."""
    mm = N
    while mm % p == 0:
        mm //= p
    return len(cyclotomic_coset_degrees(q, mm))


def is_semisimple(q, N, p):
    return N % p != 0  # gcd(N, q) == 1


# --------------------------------------------------------------------------- #
#  Best-known mu_q(n) (extension-field rank) building blocks, used by several  #
#  reference functions.  Returns (value, ref) or None.                         #
# --------------------------------------------------------------------------- #

# Cenk & Ozbudak 2010, "On multiplication in finite fields", J. Complexity 26,
# Table 1 (p. 181): upper bounds on mu_q(n) for q = 2,3,4, n = 2..18.  (This is
# the *primary* source; the same table is reproduced in the Ballet survey.)
# NB mu_3(5) is 12 here; BDEZ 2012 later improved it to 11 (handled separately).
CENK_OZBUDAK = {
    2: {
        2: 3,
        3: 6,
        4: 9,
        5: 13,
        6: 15,
        7: 22,
        8: 24,
        9: 30,
        10: 33,
        11: 39,
        12: 42,
        13: 48,
        14: 51,
        15: 54,
        16: 60,
        17: 67,
        18: 69,
    },
    3: {
        2: 3,
        3: 6,
        4: 9,
        5: 12,
        6: 15,
        7: 19,
        8: 21,
        9: 26,
        10: 27,
        11: 34,
        12: 36,
        13: 42,
        14: 45,
        15: 50,
        16: 54,
        17: 58,
        18: 62,
    },
    4: {
        2: 3,
        3: 6,
        4: 8,
        5: 11,
        6: 14,
        7: 17,
        8: 20,
        9: 23,
        10: 27,
        11: 30,
        12: 33,
        13: 37,
        14: 39,
        15: 45,
        16: 45,
        17: 53,
        18: 51,
    },
}
# Improvements over the Cenk-Ozbudak 2010 table from later exhaustive search
# (Barbulescu et al. 2012): best-known mu_q(n) upper bounds.
BDEZ_MU_UB = {(3, 5): 11}

# Short product M-hat(n) = bilinear rank of F_q[x]/(x^n), valid over any ring
# (Cenk & Ozbudak 2011, "Multiplication of polynomials modulo x^n", TCS 412,
# Table 2 diagonal, p. 3459).  Upper bounds; exact for n<=4 (cf. BDEZ).
MHAT = {1: 1, 2: 3, 3: 5, 4: 8, 5: 11, 6: 14, 7: 18, 8: 22, 9: 27, 10: 31, 11: 36}

# Full product M_q(n) = rank of (n-term) x (n-term) polynomial multiplication.
# Montgomery 2005 "five/six/seven-term Karatsuba-like formulae", valid over any
# ring, as compiled in Cenk-Ozbudak 2011 Table 1 (col. "M(n) [11]", p. 3454).
MONTGOMERY_FULL = {
    2: 3,
    3: 6,
    4: 9,
    5: 13,
    6: 17,
    7: 22,
    8: 27,
    9: 34,
    10: 39,
    11: 46,
    12: 51,
    13: 60,
    14: 66,
    15: 75,
    16: 81,
    17: 94,
    18: 102,
}
# F_2-specific full-product best-known upper bounds, as compiled in the
# Fan-Hasan 2015 survey (FFA), Table 2 "Upper bounds for M(k)" (k = terms):
# n=5,6,7 Montgomery 2005; n=8 Fan-Hasan 2007 comments; n=9 Cenk-Ozbudak 2009 /
# Oseledets 2011; n=10 FH 2007 / Oseledets.  (n<=5 proven optimal.)
F2_FULL = {2: 3, 3: 6, 4: 9, 5: 13, 6: 17, 7: 22, 8: 26, 9: 30, 10: 35}

# ---- Code-embedding lower bound data ---------------------------------------- #
# Embedding lemma (Lempel-Seroussi-Winograd 1983; Chudnovsky-Chudnovsky 1988;
# for the full product: Brockett/Brown-Dobkin, cf. Kaminski-Bshouty 1989 p.151):
# a rank-r bilinear algorithm for F_{q^n}/F_q (or for the product of two n-term
# polynomials) yields a linear [r, n, >=n]_q code, so
#       rank >= N_q(n, n) := min length of a q-ary linear code, dim n, dist n.
# N_q(n,n) >= Griesmer bound g_q(n,n) (Griesmer 1960); for q=2 the values below
# are verified LOWER bounds on N_2(k,k) from M. Grassl's tables,
# https://www.codetables.de :
# (raw pages cached under certs/codetables/, fetched by fetch_codetables.py):
#   q=2: d_max[12,5]=4 -> N=13;  d_max[14,6]=5 -> N=15;  d_max[17,7]=6 -> N=18;
#        d_max[19,8]=7 -> N=20;  d_max[25,9]=8 -> N=26;  d_max[27,10]=9 -> N=28;
#        d_max[30,11]=10, [31,11,11] exists -> N=31.
#   q=3: d_max[16,7]=6 -> N=17;  d_max[18,8]=7 -> N=19;  d_max[20,9]=8 -> N=21.
#   q=4: d_max[14,7]=6 -> N=15;  d_max[16,8]=7 -> N=17.
#   q=5: d_max[14,7]=6 -> N=15.
# (Entries equal to the Griesmer bound are omitted; code_n falls back to it.)
NQKK_EXACT = {
    (2, 5): 13,
    (2, 6): 15,
    (2, 7): 18,
    (2, 8): 20,
    (2, 9): 26,
    (2, 10): 28,
    (2, 11): 31,
    (3, 7): 17,
    (3, 8): 19,
    (3, 9): 21,
    (4, 7): 15,
    (4, 8): 17,
    (5, 7): 15,
}


def griesmer(q, k, d):
    """Griesmer 1960: any [n, k, d]_q linear code has n >= sum ceil(d/q^i)."""
    total, x = 0, d
    for i in range(k):
        total += -(-d // (q**i))
    return total


def code_n(q, k):
    """Best published lower bound on N_q(k,k), with source tag."""
    g = griesmer(q, k, k)
    if (q, k) in NQKK_EXACT and NQKK_EXACT[(q, k)] > g:
        return (NQKK_EXACT[(q, k)], "Grassl codetables.de")
    return (g, "Griesmer 1960")


def mu_ub(q, n):
    """Best published UPPER bound for mu_q(n) = rank of F_{q^n}/F_q, with source.
    Returns (value, ref) or None when nothing clean is known."""
    if n == 1:
        return (1, "trivial")
    if (q, n) in BDEZ_MU_UB:
        return (BDEZ_MU_UB[(q, n)], "Barbulescu et al. 2012, Table 3")
    if 2 * n - 1 <= q + 1:  # n <= q/2 + 1
        return (2 * n - 1, "de Groote 1983 / Toom interp.")
    if 2 * n <= q + 1 + epsilon(q):  # Shokrollahi window
        return (2 * n, "Shokrollahi 1992; Chaumine 2006")
    if q in CENK_OZBUDAK and n in CENK_OZBUDAK[q]:
        return (CENK_OZBUDAK[q][n], "Cenk-Ozbudak 2010, Table 1")
    # mu_q(n) <= M_q(n) (multiply as polynomials, reduce mod the irreducible;
    # Cenk-Ozbudak 2010 eq. (2.1)) with KB89's exact full product in its band.
    if 2 * (n - 1) > q and n - 1 <= q + 1:
        return (3 * n - 2 - q // 2, "Kaminski-Bshouty 1989 Thm 2, via mu<=M")
    return None


def mu_lb(q, n):
    """Best published LOWER bound for mu_q(n), with source (combines de Groote,
    Shokrollahi, the code-embedding bound, and specific exact values)."""
    if 2 * n - 1 <= q + 1:  # de Groote: equality regime
        return (2 * n - 1, "de Groote 1983 (Ballet survey Thm 'wdg')")
    # n > q/2 + 1  =>  strict, so mu_q(n) >= 2n.
    val, ref = 2 * n, "de Groote 1983 (>2n-1, integrality)"
    if 2 * n <= q + 1 + epsilon(q):
        ref = "Shokrollahi 1992; Chaumine 2006 (=2n)"
    cn, csrc = code_n(q, n)
    if cn > val:
        val, ref = cn, f"code embedding (LSW83/CC88) + {csrc}"
    # Stronger, specifically-proven exact lower bounds:
    if (q, n) == (2, 4):
        return (9, "Chudnovsky-Chudnovsky 1988, Ex. 3.2")
    if (q, n) == (2, 6):
        return (15, "Chudnovsky-Chudnovsky 1988, Ex. 3.3")
    if (q, n) == (3, 4) and val < 9:
        return (9, "Barbulescu et al. 2012, Table 3 (exhaustive)")
    return (val, ref)


# ---- classical CRT/interpolation upper bound for the FULL product ----------- #
# The classical algorithm (Winograd 1976; CC 1988; KB89 App. B; Fan-Hasan 2007
# comments): evaluate both t-term factors modulo pairwise-coprime prime powers
# p_i^{e_i} (including the place at infinity, treated as a linear place) of
# total degree 2t-1 and recombine by CRT.  Cost = sum of block ranks, where a
# block F_q[x]/(p^e) costs:  e=1: mu_q(deg p) (best published);  deg p=1, e>1:
# the short product (minimal 2e-1 when q >= 2e-4, else M-hat(e));  deg p>1,
# e>1: 2*deg(p)*e-1 when q >= 2*deg(p)*e-2 (AGW 1988).  Over F_2 this DP
# reproduces Fan-Hasan 2007 Table I (e.g. M_2(8) <= 26).


def count_irreducibles(q, d):
    """Number of monic irreducible polynomials of degree d over F_q (Gauss)."""

    def moebius(x):
        out, i = 1, 2
        while i * i <= x:
            if x % i == 0:
                x //= i
                if x % i == 0:
                    return 0
                out = -out
            i += 1
        return -out if x > 1 else out

    total, e = 0, 1
    while e <= d:
        if d % e == 0:
            total += moebius(e) * q ** (d // e)
        e += 1
    return total // d


def i_q_max_factors(q, d):
    """i_q(d): the maximal number of distinct irreducible factors of a degree-d
    polynomial over F_q (greedily take the smallest-degree irreducibles).
    Kaminski-Bshouty 1989, Appendix A."""
    count, used, deg = 0, 0, 1
    while deg <= d:
        can = min(count_irreducibles(q, deg), (d - used) // deg)
        if can <= 0:
            break
        count += can
        used += can * deg
        deg += 1
    return count


def _crt_block_cost(q, d, e):
    """Closed-form block costs available to the classical CRT algorithm
    (deliberately recursion-free)."""
    if e == 1:
        if d == 1:
            return 1
        b = mu_ub(q, d)
        return b[0] if b else None
    if d == 1:
        if q >= 2 * e - 4:
            return 2 * e - 1
        return MHAT.get(e)
    return (2 * d * e - 1) if q >= 2 * d * e - 2 else None


def crt_interp_full_ub(q, t):
    """Min cost of the classical CRT/interpolation algorithm for the full
    product of two t-term polynomials over F_q (valid upper bound), or None."""
    D = 2 * t - 1
    INF = 10**9
    best = [INF] * (D + 1)  # best[w] = min cost to cover total degree w
    best[0] = 0
    for d in range(1, D + 1):
        supply = (q + 1) if d == 1 else count_irreducibles(q, d)
        supply = min(supply, D // d)
        if supply <= 0:
            continue
        # options for ONE prime of degree d: multiplicity e -> (weight, cost)
        opts = []
        for e in range(1, D // d + 1):
            cst = _crt_block_cost(q, d, e)
            if cst is not None:
                opts.append((d * e, cst))
        cur = best[:]
        for _ in range(supply):  # each prime used at most once
            nxt = cur[:]
            for w in range(D + 1):
                if cur[w] == INF:
                    continue
                for wt, cst in opts:
                    w2 = min(D, w + wt)
                    if cur[w] + cst < nxt[w2]:
                        nxt[w2] = cur[w] + cst
            cur = nxt
        best = cur
    return best[D] if best[D] < INF else None


def full_ub(q, t):
    """Best published upper bound for the t-term full product M_q(t), with
    source."""
    cands = []
    if 2 * t - 1 <= q + 1:
        cands.append((2 * t - 1, "Toom 1963 / interpolation"))
    if 2 * (t - 1) > q and t - 1 <= q + 1:
        cands.append((3 * t - 2 - q // 2, "Kaminski-Bshouty 1989, Thm 2"))
    if q == 2 and t in F2_FULL:
        cands.append((F2_FULL[t], "Fan-Hasan 2015 survey, Table 2"))
    if t in MONTGOMERY_FULL:
        cands.append((MONTGOMERY_FULL[t], "Montgomery 2005"))
    v = crt_interp_full_ub(q, t)
    if v is not None:
        cands.append((v, "CRT/interpolation (Winograd 1976; cf. FH 2007)"))
    return min(cands) if cands else None


# ---- per-block bounds for F_q[x]/(p(x)^e), deg p = d  (used for 'cyclic') --- #
# Minimal-rank criteria over F_q:
#  * e = 1, d = 1:  F_q, rank 1.
#  * e = 1, d >= 2: field F_{q^d}: rank 2d-1 iff q >= 2d-2 (Winograd/de Groote).
#  * d = 1, e >= 2: F_q[y]/(y^e): rank 2e-1 iff q >= 2e-4
#       (Averbuch-Galil-Winograd Part II, TCS 86 (1991); Blaser 2005, Thm 22).
#  * d >= 2, e >= 2: F_q[x]/(p^e): rank 2de-1 iff q >= 2de-2
#       (AGW Part I, TCS 58 (1988), Thms 5.1/5.1'; Blaser 2005, Thm 35).
# Failure of the criterion implies rank >= 2*dim (strict + integrality).


def block_lb(q, d, e):
    """Best published lower bound for rank of F_q[x]/(p^e), deg p = d."""
    dim = d * e
    if dim == 1:
        return 1
    if e == 1:
        return mu_lb(q, d)[0]
    if d == 1:
        if q >= 2 * e - 4:
            return 2 * e - 1  # minimal (AGW II / Blaser)
        if (q, e) == (2, 4):
            return 8  # BDEZ 2012, Table 4 (exact)
        if (q, e) == (2, 5):
            return 10  # BDEZ 2012 (rank > 9)
        return 2 * e  # not minimal
    return (2 * dim - 1) if q >= 2 * dim - 2 else 2 * dim


def block_ub(q, d, e):
    """Best published upper bound for rank of F_q[x]/(p^e), with source, or
    None.  Always also a candidate: FOLD the full product of two (de)-term
    polynomials and reduce mod p^e (linear), so rank <= M_q(de)."""
    dim = d * e
    if dim == 1:
        return (1, "trivial")
    fold = full_ub(q, dim)
    fold = (fold[0], "fold of M_q(%d) [%s]" % (dim, fold[1])) if fold else None
    cands = [fold] if fold else []
    if e == 1:
        b = mu_ub(q, d)
        if b:
            cands.append(b)
    elif d == 1:
        if q >= 2 * e - 4:
            cands.append((2 * e - 1, "AGW 1991 / Blaser 2005 Thm 22 (minimal)"))
        elif e in MHAT:
            cands.append((MHAT[e], "Cenk-Ozbudak 2011, Table 2"))
    else:
        if q >= 2 * dim - 2:
            cands.append((2 * dim - 1, "AGW 1988, Thm 5.1' (minimal)"))
    return min(cands) if cands else None


def cyclic_blocks(q, N, p):
    """Distinct prime-power blocks of x^N - 1 over F_q: list of (d, e)."""
    mm, a = N, 0
    while mm % p == 0:
        mm //= p
        a += 1
    e = p**a
    return [(d, e) for d in cyclotomic_coset_degrees(q, mm)]


# --------------------------------------------------------------------------- #
#  REFERENCE FUNCTIONS.  Each returns a list of Claim for the given cell.      #
#  cell = dict(problem, p, m, q, n).                                           #
# --------------------------------------------------------------------------- #


def ref_fiduccia_zalcstein_1971(c):
    """Fiduccia & Zalcstein 1971 (cf. Buergisser-Clausen-Shokrollahi Prop.14.47;
    Ballet survey Thm 'old'): computing the 2n-1 coefficients of a product of two
    degree-(n-1) polynomials needs >= 2n-1 multiplications.  Lower bound for the
    'full' product and (mod an irreducible) for 'extfield'."""
    n = c["n"]
    if c["problem"] == "full":
        return [Claim("lb", 2 * n - 1, "Fiduccia-Zalcstein 1971", 1)]
    if c["problem"] == "extfield":
        return [Claim("lb", 2 * n - 1, "Fiduccia-Zalcstein 1971", 1)]
    return []


def ref_winograd_1976_cor1(c):
    """S. Winograd 1976, 'Some bilinear forms whose multiplicative complexity
    depends on the field of constants', Cor. 1 / Thm 5: multiplication modulo a
    polynomial P with k distinct irreducible factors needs >= 2*deg(P) - k
    multiplications.  For 'cyclic', P = x^N - 1."""
    if c["problem"] != "cyclic":
        return []
    q, N, p = c["q"], c["n"], c["p"]
    k = num_distinct_irreducible_factors(q, N, p)
    return [Claim("lb", 2 * N - k, "Winograd 1976, Cor. 1", 1)]


def ref_degroote_1983(c):
    """H. de Groote 1983 + S. Winograd 1979 (Ballet survey Thm 'wdg'):
    mu_q(n) >= 2n-1 with equality iff n <= q/2 + 1.  Gives the exact value
    2n-1 in that range (lb and ub), and the strict lower bound 2n beyond it."""
    if c["problem"] != "extfield":
        return []
    q, n = c["q"], c["n"]
    out = []
    if 2 * n - 1 <= q + 1:  # exact 2n-1
        out.append(Claim("lb", 2 * n - 1, "de Groote 1983 (=2n-1, n<=q/2+1)", 3))
        out.append(Claim("ub", 2 * n - 1, "de Groote 1983 / Toom interp.", 3))
    else:  # strict: >= 2n
        out.append(Claim("lb", 2 * n, "de Groote 1983 (>2n-1, integrality)", 2))
    return out


def ref_shokrollahi_1992(c):
    """Shokrollahi 1992 (elliptic curves) + Chaumine: mu_q(n) = 2n for
    q/2+1 < n <= (q+1+eps(q))/2.  Exact value (lb and ub)."""
    if c["problem"] != "extfield":
        return []
    q, n = c["q"], c["n"]
    if (2 * n - 1 > q + 1) and (2 * n <= q + 1 + epsilon(q)):
        return [
            Claim("lb", 2 * n, "Shokrollahi 1992; Chaumine 2006", 3),
            Claim("ub", 2 * n, "Shokrollahi 1992; Chaumine 2006", 3),
        ]
    return []


def ref_chudnovsky_1988(c):
    """D.V. & G.V. Chudnovsky 1988, 'Algebraic complexities and algebraic curves
    over finite fields', Examples 3.2/3.3: mu_2(4)=9 and mu_2(6)=15 (exact)."""
    if c["problem"] != "extfield":
        return []
    q, n = c["q"], c["n"]
    exact = {(2, 4): 9, (2, 6): 15}
    if (q, n) in exact:
        v = exact[(q, n)]
        return [
            Claim("lb", v, "Chudnovsky-Chudnovsky 1988, Ex. 3.2/3.3", 4),
            Claim("ub", v, "Chudnovsky-Chudnovsky 1988, Ex. 3.2/3.3", 4),
        ]
    return []


def ref_toom_interpolation_full(c):
    """A. Toom 1963 / S. Winograd 1979 evaluation-interpolation: a product of two
    degree-(n-1) polynomials over F_q can be computed in 2n-1 multiplications
    whenever there are enough points, i.e. 2n-1 <= q+1 (n <= q/2+1).  Matching
    Fiduccia-Zalcstein, this is the exact value M_q(n)=2n-1 in that range."""
    if c["problem"] != "full":
        return []
    q, n = c["q"], c["n"]
    if 2 * n - 1 <= q + 1:
        return [
            Claim("ub", 2 * n - 1, "Toom 1963 / interpolation", 3),
            Claim("lb", 2 * n - 1, "Fiduccia-Zalcstein 1971 (exact, n<=q/2+1)", 3),
        ]
    return []


def ref_cenk_ozbudak_2010(c):
    """M. Cenk & F. Ozbudak 2010, "On multiplication in finite fields",
    J. Complexity 26, Table 1 (p. 181): best-known upper bounds on mu_q(n) for
    q in {2,3,4}.  Used where de Groote/Shokrollahi do not already give an exact
    value.  (mu_3(5)=12 here is superseded by BDEZ 2012's 11.)"""
    if c["problem"] != "extfield":
        return []
    q, n = c["q"], c["n"]
    if 2 * n - 1 <= q + 1 or 2 * n <= q + 1 + epsilon(q):
        return []
    if q in CENK_OZBUDAK and n in CENK_OZBUDAK[q]:
        return [Claim("ub", CENK_OZBUDAK[q][n], "Cenk-Ozbudak 2010, Table 1", 2)]
    return []


def ref_montgomery_full(c):
    """P. Montgomery 2005, "Five, six, and seven-term Karatsuba-like formulae"
    (IEEE TC 54), as compiled in Cenk-Ozbudak 2011 Table 1.  Karatsuba-like full
    products valid over ANY ring, hence an upper bound on M_q(n) for every q."""
    if c["problem"] != "full":
        return []
    n = c["n"]
    if n in MONTGOMERY_FULL:
        return [
            Claim(
                "ub",
                MONTGOMERY_FULL[n],
                "Montgomery 2005 (via Cenk-Ozbudak 2011, Table 1)",
                2,
            )
        ]
    return []


def ref_f2_full(c):
    """Best-known F_2 full-product upper bounds as compiled in the Fan-Hasan
    2015 survey (FFA), Table 2 'Upper bounds for M(k)': n=5,6,7 Montgomery 2005
    (13/17/22); n=8 Fan-Hasan 2007 comments (26, below Montgomery's 27);
    n=9 Cenk-Ozbudak 2009 / Oseledets 2011 (30); n=10 (35)."""
    if c["problem"] != "full" or c["q"] != 2:
        return []
    n = c["n"]
    if n in F2_FULL:
        return [Claim("ub", F2_FULL[n], "Fan-Hasan 2015 survey, Table 2", 3)]
    return []


# --- Barbulescu-Detrey-Estibals-Zimmermann 2012 exhaustive search (F_2, F_3) -- #
# Tables 1-4 of "Finding Optimal Formulae for Bilinear Maps".  A non-symmetric
# row that found no solution at k proves rank > k (rigorous lower bound rank>=k+1
# over F_2/F_3); a row with solutions at the smallest tried k proves the exact
# rank; a symmetric "(Sym.)" row with a solution gives a valid upper bound.
BDEZ_EXACT = {  # (problem, q, n): rank  (lb == ub)
    ("full", 2, 2): 3,
    ("full", 2, 3): 6,
    ("full", 2, 4): 9,
    ("full", 2, 5): 13,
    ("full", 3, 2): 3,
    ("full", 3, 3): 6,
    ("full", 3, 4): 9,
    ("full", 3, 5): 12,
    ("extfield", 2, 2): 3,
    ("extfield", 2, 3): 6,
    ("extfield", 2, 4): 9,
    ("extfield", 3, 2): 3,
    ("extfield", 3, 3): 6,
    ("extfield", 3, 4): 9,
    ("cyclic", 2, 2): 3,
    ("cyclic", 2, 3): 4,
    ("cyclic", 2, 4): 8,
    ("cyclic", 2, 5): 10,
    ("cyclic", 3, 2): 2,
    ("cyclic", 3, 3): 5,
    ("cyclic", 3, 4): 5,
}
BDEZ_LB = {  # (problem, q, n): proven rank >= value (no soln at value-1)
    ("full", 2, 6): 15,  # 6x6 over F2: no formula with 14 -> rank > 14
}
BDEZ_UB = {  # (problem, q, n): best symmetric formula found
    ("full", 2, 6): 17,
    ("full", 2, 7): 22,
    ("extfield", 2, 5): 13,
    ("extfield", 2, 6): 15,
    ("extfield", 3, 5): 11,
    ("cyclic", 2, 6): 12,
    ("cyclic", 2, 7): 13,
    ("cyclic", 3, 5): 10,
    ("cyclic", 3, 6): 10,
}


def ref_bdez_2012(c):
    key = (c["problem"], c["q"], c["n"])
    out = []
    if key in BDEZ_EXACT:
        v = BDEZ_EXACT[key]
        out.append(Claim("lb", v, "Barbulescu et al. 2012 (exhaustive, exact)", 4))
        out.append(Claim("ub", v, "Barbulescu et al. 2012 (exhaustive, exact)", 4))
    if key in BDEZ_LB:
        out.append(Claim("lb", BDEZ_LB[key], "Barbulescu et al. 2012 (exhaustive)", 4))
    if key in BDEZ_UB:
        out.append(
            Claim("ub", BDEZ_UB[key], "Barbulescu et al. 2012 (symmetric formula)", 4)
        )
    return out


def ref_kb89_general_lb(c):
    """Kaminski-Bshouty 1989 (J. ACM 36), Theorem 1 (proof, p. 158): for q >= 3,
    M_q(deg) >= 3*deg + 1 - i_q(deg), where i_q(d) is the maximal number of
    distinct irreducible factors of a degree-d polynomial over F_q.  (deg = n-1
    in our n-term convention; for q=2 the paper instead uses the coding bound,
    handled by ref_code_embedding.)"""
    if c["problem"] != "full" or c["q"] < 3:
        return []
    d = c["n"] - 1
    if d < 1:
        return []
    v = 3 * d + 1 - i_q_max_factors(c["q"], d)
    return [Claim("lb", v, "Kaminski-Bshouty 1989, Thm 1 (3d+1-i_q(d))", 4)]


def ref_mu_le_full_product(c):
    """mu_q(n) <= M_q(n): to multiply in F_{q^n}=F_q[x]/(P) one may multiply the
    two representatives as n-term polynomials and reduce modulo P (a linear,
    cost-free step).  Cenk-Ozbudak 2010 eq. (2.1).  Combined with the best
    published M_q(n) (here the classical CRT/interpolation construction), this
    upper-bounds the tensor rank of every extension -- in particular for q and n
    outside the de Groote/Shokrollahi/Cenk-Ozbudak ranges."""
    if c["problem"] != "extfield":
        return []
    b = full_ub(c["q"], c["n"])
    if b is None:
        return []
    return [Claim("ub", b[0], "mu(n)<=M_q(n); " + b[1], 1)]


def ref_crt_interp_full(c):
    """Classical CRT / evaluation-interpolation algorithm for the full product
    of two n-term polynomials over F_q (Winograd 1976; Chudnovsky-Chudnovsky
    1988; KB89 App. B; the construction tabulated by Fan-Hasan 2007).  Upper
    bound = crt_interp_full_ub(q, n); valid for every q and n."""
    if c["problem"] != "full":
        return []
    v = crt_interp_full_ub(c["q"], c["n"])
    if v is None:
        return []
    return [Claim("ub", v, "CRT/interpolation (Winograd 1976; cf. FH 2007)", 2)]


def ref_kaminski_bshouty_1989(c):
    """Kaminski & Bshouty 1989 (J. ACM 36), Theorem 2: for any q and
    q/2 < deg <= q+1 (deg = n-1 in our n-term convention), the full product has
    EXACTLY  3*deg + 1 - floor(q/2)  multiplications, i.e. 3n - 2 - floor(q/2).
    (For deg <= q/2 the classical value 2n-1 holds, KB89 p. 151.)
    For 'extfield', mu_q(n) <= M_q(n) (multiply as polynomials, then reduce
    modulo the irreducible -- linear; Cenk-Ozbudak 2010 eq. (2.1)) turns the
    same exact value into an upper bound."""
    if c["problem"] not in ("full", "extfield"):
        return []
    q, n = c["q"], c["n"]
    deg = n - 1
    if not (2 * deg > q and deg <= q + 1):
        return []
    v = 3 * n - 2 - q // 2
    if c["problem"] == "full":
        return [
            Claim("lb", v, "Kaminski-Bshouty 1989, Thm 2 (exact)", 5),
            Claim("ub", v, "Kaminski-Bshouty 1989, Thm 2 (exact)", 5),
        ]
    return [Claim("ub", v, "Kaminski-Bshouty 1989 Thm 2, via mu<=M", 1)]


def ref_code_embedding(c):
    """Code-embedding lower bound: a rank-r algorithm for F_{q^n}/F_q (LSW 1983;
    CC 1988) or for the n-term full product (Brockett/Brown-Dobkin, cf. KB89
    p. 151) yields an [r, n, >=n]_q linear code, so rank >= N_q(n,n).
    N values: Griesmer 1960 bound, sharpened by Grassl's codetables.de for q=2
    (e.g. no [12,5,5]_2 (Fontaine-Peterson 1959) -> mu_2(5) >= 13)."""
    if c["problem"] not in ("extfield", "full"):
        return []
    q, n = c["q"], c["n"]
    if n < 2:
        return []
    v, src = code_n(q, n)
    return [Claim("lb", v, f"code embedding (LSW83/CC88) + {src}", 2)]


def ref_blaser_cyclic_lb(c):
    """Blaser 2005 (SIAM JC 34), proof of Lemma 26 / Lemma 11: for a product of
    local algebras, R(prod B_i) >= sum_{i != j} (2 dim B_i - 1) + R(B_j).
    Applied to F_q[x]/(x^N-1) = prod F_q[x]/(p_i^e) with the best published
    lower bound for one block B_j (AGW 1988/1991 minimal-rank criteria, the
    code-embedding bound, or exact values), maximised over j.  Emitted only
    when it beats Winograd's 2N - k."""
    if c["problem"] != "cyclic":
        return []
    q, N, p = c["q"], c["n"], c["p"]
    blocks = cyclic_blocks(q, N, p)
    k = len(blocks)
    total_as = sum(2 * d * e - 1 for d, e in blocks)  # = 2N - k
    best = max(total_as - (2 * d * e - 1) + block_lb(q, d, e) for d, e in blocks)
    if best > 2 * N - k:
        return [Claim("lb", best, "Blaser 2005, Lem. 26 + AGW/code block bounds", 3)]
    return []


def ref_wagh_morgera_1983(c):
    """Wagh & Morgera 1983 (IEEE IT-29), Tables I-IV: constructed cyclic-
    convolution algorithms over GF(2) and GF(3) (primary source predating
    Morgera 1990); explicitly tabulated small lengths."""
    if c["problem"] != "cyclic":
        return []
    WM = {
        (2, 3): 4,
        (2, 5): 10,
        (2, 7): 13,
        (2, 9): 22,
        (3, 2): 2,
        (3, 4): 5,
        (3, 5): 10,
        (3, 7): 19,
        (3, 8): 11,
        (3, 10): 20,
    }
    key = (c["q"], c["n"])
    if key in WM:
        return [Claim("ub", WM[key], "Wagh-Morgera 1983, Tables I-IV", 5)]
    return []


def ref_crt_cyclic(c):
    """CRT decomposition of F_q[x]/(x^N-1) = prod_i F_q[x]/(p_i^e), the cyclic-
    convolution method of Morgera 1990 (Eq. 3) [for prime q; Winograd 1976 /
    LSW 1983 in general].  Upper bound = sum of best published block ranks
    (block_ub: de Groote/Shokrollahi/Cenk-Ozbudak 2010 for field blocks,
    AGW 1988/1991 + Blaser 2005 minimal-rank values or Cenk-Ozbudak 2011 M-hat
    for repeated factors); emitted only when every block has a published value."""
    if c["problem"] != "cyclic":
        return []
    q, N, p = c["q"], c["n"], c["p"]
    total, srcs = 0, set()
    for d, e in cyclic_blocks(q, N, p):
        b = block_ub(q, d, e)
        if b is None:
            return []
        total += b[0]
        srcs.add(b[1].split("(")[0].split(",")[0].strip())
    method = "Morgera 1990" if c["m"] == 1 else "Winograd 1976/LSW 1983"
    label = f"cyclic CRT ({method}); " + "; ".join(sorted(srcs))
    return [Claim("ub", total, label, 2)]


REFERENCE_FUNCTIONS = [
    ref_fiduccia_zalcstein_1971,
    ref_winograd_1976_cor1,
    ref_degroote_1983,
    ref_shokrollahi_1992,
    ref_chudnovsky_1988,
    ref_toom_interpolation_full,
    ref_crt_interp_full,
    ref_mu_le_full_product,
    ref_kaminski_bshouty_1989,
    ref_kb89_general_lb,
    ref_code_embedding,
    ref_cenk_ozbudak_2010,
    ref_montgomery_full,
    ref_f2_full,
    ref_bdez_2012,
    ref_blaser_cyclic_lb,
    ref_wagh_morgera_1983,
    ref_crt_cyclic,
]

# --------------------------------------------------------------------------- #
#  Combine                                                                     #
# --------------------------------------------------------------------------- #


def best_bounds(cell):
    claims = []
    for f in REFERENCE_FUNCTIONS:
        claims.extend(f(cell))
    lbs = [c for c in claims if c.kind == "lb"]
    ubs = [c for c in claims if c.kind == "ub"]
    lb = ub = None
    if lbs:
        m = max(c.value for c in lbs)
        lb = max((c for c in lbs if c.value == m), key=lambda c: c.prio)
    if ubs:
        m = min(c.value for c in ubs)
        ub = max((c for c in ubs if c.value == m), key=lambda c: c.prio)
    return lb, ub


def main():
    rows = []
    for problem in PROBLEMS:
        for q, p, m in prime_powers(Q_MAX):
            for n in range(1, N_MAX + 1):
                cell = dict(problem=problem, p=p, m=m, q=q, n=n)
                lb, ub = best_bounds(cell)
                rows.append(
                    {
                        "problem": problem,
                        "p": p,
                        "m": m,
                        "q": q,
                        "n": n,
                        "lb": lb.value if lb else "-",
                        "lb_ref": lb.ref if lb else "-",
                        "ub": ub.value if ub else "-",
                        "ub_ref": ub.ref if ub else "-",
                    }
                )
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["problem", "p", "m", "q", "n", "lb", "lb_ref", "ub", "ub_ref"],
        )
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()


# --------------------------------------------------------------------------- #
#  REFERENCES (one paper per citation label used above)                        #
# --------------------------------------------------------------------------- #
# Fiduccia-Zalcstein 1971  C.M. Fiduccia, Y. Zalcstein, "Algebras having linear
#     multiplicative complexities", JACM 24 (1977) 311-331 (lower bound 2n-1;
#     also Buergisser-Clausen-Shokrollahi, "Algebraic Complexity Theory",
#     Prop. 14.47).
# Toom 1963  A.L. Toom, "The complexity of a scheme of functional elements...",
#     Soviet Math. Dokl. 3 (1963) 714-716 (interpolation upper bound 2n-1).
# Winograd 1976  S. Winograd, "Some bilinear forms whose multiplicative
#     complexity depends on the field of constants", Math. Systems Theory 10
#     (1977) 169-180 (Cor. 1 / Thm 5: mult mod P needs >= 2 deg(P) - k).
# Winograd 1979  S. Winograd, "On multiplication in algebraic extension fields",
#     Theoret. Comput. Sci. 8 (1979) 359-377.
# de Groote 1983  H.F. de Groote, "Characterization of division algebras of
#     minimal rank...", SIAM J. Comput. 12 (1983) 101-117 (mu_q(n)=2n-1 iff
#     n <= q/2+1).
# Lempel-Seroussi-Winograd 1983  TCS 22 (1983) 285-296.
# Chudnovsky-Chudnovsky 1988  D.V./G.V. Chudnovsky, "Algebraic complexities and
#     algebraic curves over finite fields", J. Complexity 4 (1988) 285-316
#     (Ex. 3.2/3.3: mu_2(4)=9, mu_2(6)=15).
# Morgera 1990  S.D. Morgera, "Multiplicative complexity of bilinear algorithms
#     for cyclic convolution over finite fields", Multidim. Syst. Signal Process.
#     1 (1990) 99-111 (Eq. 3 / Thm 1; Tables 1-2: cyclic-convolution CRT method
#     over GF(2), GF(3)).
# Shokrollahi 1992  M.A. Shokrollahi, "Optimal algorithms for multiplication in
#     certain finite fields using elliptic curves", SIAM J. Comput. 21 (1992)
#     1193-1198; Chaumine 2006 (AGCT) (mu_q(n)=2n in the elliptic window).
# Montgomery 2005  P.L. Montgomery, "Five, six, and seven-term Karatsuba-like
#     formulae", IEEE Trans. Comput. 54 (2005) 362-369 (full products over any
#     ring; values via Cenk-Ozbudak 2011 Table 1).
# Cenk-Ozbudak 2009  M. Cenk, F. Ozbudak, "Improved polynomial multiplication
#     formulas over F_2 using CRT", IEEE Trans. Comput. 58 (2009) (F_2 full
#     products; values via CO 2011 Table 1).  Fan-Hasan: H. Fan, M.A. Hasan.
# Cenk-Ozbudak 2010  M. Cenk, F. Ozbudak, "On multiplication in finite fields",
#     J. Complexity 26 (2010) 172-186 (Table 1: mu_q(n), q=2,3,4, n=2..18).
# Cenk-Ozbudak 2011  M. Cenk, F. Ozbudak, "Multiplication of polynomials modulo
#     x^n", Theoret. Comput. Sci. 412 (2011) 3451-3462 (Table 1: full M_q(n);
#     Table 2: short product M-hat(n) = rank of F_q[x]/(x^n)).
# Barbulescu-Detrey-Estibals-Zimmermann 2012  "Finding Optimal Formulae for
#     Bilinear Maps", WAIFI 2012, LNCS 7369, 168-186 (exhaustive/symmetric
#     search over F_2, F_3: Tables 1-4).
# Wagh-Morgera 1983  "A new structured design method for convolutions over
#     finite fields", IEEE Trans. Inform. Theory IT-29 (1983) 583-595
#     (Tables I-IV: cyclic convolution over GF(2)/GF(3)).
# Averbuch-Galil-Winograd 1988/1991  "Classification of all the minimal bilinear
#     algorithms for ... mod Q(u)^l", Part I: TCS 58 (1988) 17-56 (deg Q >= 2:
#     minimal 2n-1 iff q >= 2n-2, Thms 5.1/5.1'); Part II: TCS 86 (1991) 143-203
#     (Q = u: F_q[u]/(u^n) minimal iff q >= 2n-4, cf. Blaser 2005 Thm 22).
# Kaminski-Bshouty 1989  "Multiplicative complexity of polynomial multiplication
#     over finite fields", J. ACM 36 (1989) 150-170 (Thm 2: exact
#     M = 3deg+1-floor(q/2) for q/2 < deg <= q+1; coding lower bound p. 151).
# Blaser 2005  "A complete characterization of the algebras of minimal bilinear
#     complexity", SIAM J. Comput. 34 (2005) 277-298 (Thm 35, Thm 22, Lem. 26:
#     per-factor minimal-rank criteria; product minimal iff every factor is).
# Griesmer 1960  "A bound for error-correcting codes", IBM J. Res. Dev. 4
#     (1960) 532-542 (length bound used in the code-embedding lower bound).
# Grassl codetables  M. Grassl, "Bounds on the minimum distance of linear
#     codes", https://www.codetables.de (exact N_2(k,k) for k = 5..8; entries
#     trace to Fontaine-Peterson 1959 a.o.).
# Fan-Hasan 2015  "A survey of some recent bit-parallel GF(2^n) multipliers",
#     Finite Fields Appl. 32 (2015) 5-43 (Table 2: best-known M(k) over F_2).
# Fan-Hasan 2007  "Comments on 'Five, six, and seven-term Karatsuba-like
#     formulae'", IEEE Trans. Comput. 56 (2007) 716-717 (M_2(8) <= 26).
# Oseledets 2011  "Improved n-term Karatsuba-like formulas in GF(2)", IEEE
#     Trans. Comput. 60 (2011) 1212-1216.
# (Checked, no usable small-parameter values for these cells: Ballet-Pieltant
#     2011 (J. Compl. 27, asymptotic only); Randriambololona 2012 (J. Compl. 28,
#     algebra bounds mu_q(m,l) only); Ballet-Bonnecaze-Pacifico 2021
#     (arXiv:2007.16082, projective-line Chudnovsky: values match or exceed the
#     tables above; mu_4(3)=5 already de Groote); Alexeev-Forbes-Tsimerman 2011;
#     Lavrauw-Sheekey 2021 (arXiv:2102.01997, semifield ranks: confirms
#     mu_2(4)=9, mu_3(4)=9); Bshouty-Kaminski 1990 (straight-line model).
