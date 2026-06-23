#!/usr/bin/env python3
import math
import csv

# ==========================================
# Algebraic Core Logic
# ==========================================


def get_factors_info(n, q, p, is_negacyclic=False):
    """
    Computes `k` (number of distinct irreducible factors) and
    `max_d` (maximum dimension among local algebras) for x^n ± 1 over F_q.
    This works by evaluating orbits under multiplication by q.
    """
    v = 0
    n0 = n
    # Factor out the characteristic p
    while n0 % p == 0:
        v += 1
        n0 //= p

    e = p**v  # Multiplicity of each root

    # In characteristic 2, x^n + 1 == x^n - 1
    if is_negacyclic and p == 2:
        is_negacyclic = False

    visited = set()
    k = 0
    max_d = 0

    if not is_negacyclic:
        # Roots of x^n0 - 1 are z^a for a in 0..n0-1
        for a in range(n0):
            if a not in visited:
                k += 1
                orbit_size = 0
                curr = a
                while curr not in visited:
                    visited.add(curr)
                    orbit_size += 1
                    curr = (curr * q) % n0
                D_i = e * orbit_size
                if D_i > max_d:
                    max_d = D_i
    else:
        # Roots of x^n0 + 1 are z^a for odd a in 1, 3, ..., 2n0-1
        for a in range(1, 2 * n0, 2):
            if a not in visited:
                k += 1
                orbit_size = 0
                curr = a
                while curr not in visited:
                    visited.add(curr)
                    orbit_size += 1
                    curr = (curr * q) % (2 * n0)
                D_i = e * orbit_size
                if D_i > max_d:
                    max_d = D_i

    return k, max_d


# ==========================================
# Lemma / Theorem Functions
# ==========================================


def winograd_lb(n):
    return 2 * n - 1, "Winograd (1977)"


def alder_strassen_lb(n, k):
    return 2 * n - k, "Alder-Strassen (1981)"


def winograd_tfull_ub(n, q):
    if q >= 2 * n - 2:
        return 2 * n - 1, "Winograd (1977)"
    return float("inf"), ""


def kaminski_bshouty_tfull(n, q):
    """Exact bound for T_full when n-2 <= q < 2n-2."""
    if n - 2 <= q < 2 * n - 2:
        val = 3 * n - 2 - (q // 2)
        return val, "Kaminski & Bshouty (1989)"
    return float("inf"), ""


def lempel_winograd_ub(n, q, max_d, k):
    """Exact Lempel-Winograd bound for generic quotient rings if q is large enough."""
    if q >= 2 * max_d - 2:
        return 2 * n - k, "Lempel-Winograd (1983)"
    return float("inf"), ""


def ballet_text_ub(n, q):
    if q >= 2 * n - 2:
        return 2 * n - 1, "Ballet et al. (2021)"
    return float("inf"), ""


# ==========================================
# Main Aggregation Logic
# ==========================================


def calculate_bounds():
    results = []
    # Primes p and degrees m for q = p^m <= 16
    domain = [
        (2, 1),
        (2, 2),
        (2, 3),
        (2, 4),
        (3, 1),
        (3, 2),
        (5, 1),
        (7, 1),
        (11, 1),
        (13, 1),
    ]

    def add_result(prob, lb, lb_ref, ub, ub_ref):
        results.append(
            {
                "problem": prob,
                "p": p,
                "m": m,
                "q": q,
                "n": n,
                "lb": lb,
                "lb_ref": lb_ref,
                "ub": ub if ub != float("inf") else "",
                "ub_ref": ub_ref,
            }
        )

    for p, m in domain:
        q = p**m
        for n in range(1, 12):
            # --- 1. T_full ---
            lb_full, lb_ref_full = winograd_lb(n)
            ub_full, ub_ref_full = float("inf"), ""

            w_ub, w_ref = winograd_tfull_ub(n, q)
            if w_ub < ub_full:
                ub_full, ub_ref_full = w_ub, w_ref

            kb_val, kb_ref = kaminski_bshouty_tfull(n, q)
            if kb_val != float("inf"):
                lb_full, lb_ref_full = kb_val, kb_ref
                ub_full, ub_ref_full = kb_val, kb_ref

            add_result("full", lb_full, lb_ref_full, ub_full, ub_ref_full)

            # --- 3. T_cyc ---
            k_cyc, max_d_cyc = get_factors_info(n, q, p, is_negacyclic=False)
            lb_cyc, lb_ref_cyc = alder_strassen_lb(n, k_cyc)
            ub_cyc, ub_ref_cyc = ub_full, ub_ref_full

            lw_val_cyc, lw_ref_cyc = lempel_winograd_ub(n, q, max_d_cyc, k_cyc)
            if lw_val_cyc != float("inf"):
                ub_cyc, ub_ref_cyc = lw_val_cyc, lw_ref_cyc

            add_result("cyclic", lb_cyc, lb_ref_cyc, ub_cyc, ub_ref_cyc)

            # --- 4. T_negacyc ---
            k_neg, max_d_neg = get_factors_info(n, q, p, is_negacyclic=True)
            lb_neg, lb_ref_neg = alder_strassen_lb(n, k_neg)
            ub_neg, ub_ref_neg = ub_full, ub_ref_full

            lw_val_neg, lw_ref_neg = lempel_winograd_ub(n, q, max_d_neg, k_neg)
            if lw_val_neg != float("inf"):
                ub_neg, ub_ref_neg = lw_val_neg, lw_ref_neg

            add_result("negacyclic", lb_neg, lb_ref_neg, ub_neg, ub_ref_neg)

            # --- 2. T_trunc ---
            # Algebra F_q[x]/(x^n) has exactly 1 maximal ideal -> k=1
            lb_trunc, lb_ref_trunc = alder_strassen_lb(n, k=1)
            ub_trunc, ub_ref_trunc = ub_full, ub_ref_full  # Caps at T_full

            lw_val_trunc, lw_ref_trunc = lempel_winograd_ub(n, q, max_d=n, k=1)
            if lw_val_trunc < ub_trunc:
                ub_trunc, ub_ref_trunc = lw_val_trunc, lw_ref_trunc

            add_result("truncated", lb_trunc, lb_ref_trunc, ub_trunc, ub_ref_trunc)

            # --- 5. T_ext ---
            lb_ext, lb_ref_ext = winograd_lb(n)
            ub_ext, ub_ref_ext = ub_full, ub_ref_full

            ballet_val, ballet_ref = ballet_text_ub(n, q)
            if ballet_val != float("inf"):
                ub_ext, ub_ref_ext = ballet_val, ballet_ref

            add_result("extfield", lb_ext, lb_ref_ext, ub_ext, ub_ref_ext)

    return results


def write_csv(results, filename="prev_results_gemini.csv"):
    keys = ["problem", "p", "m", "q", "n", "lb", "lb_ref", "ub", "ub_ref"]
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    res = calculate_bounds()
    write_csv(res)
    print(f"Successfully generated prev_results_gemini.csv with {len(res)} records.")
