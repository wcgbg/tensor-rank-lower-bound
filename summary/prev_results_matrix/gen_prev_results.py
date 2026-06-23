#!/usr/bin/env python3
"""Generate prev_results.csv: previously-published best-known lower bounds
for the bilinear complexity (tensor rank) of <n0,n1,n2> matrix multiplication over
F_2, the matrix analogue of the polynomial prev_results_polynomial/gen_prev_results.py.

Each reference is one Python function that, given (n0,n1,n2), returns
(value, -year, tag) -- the negated year makes max() break ties toward the
earliest-published reference.  rank_lower_bound() maximises over all four
references and all six permutations of the three sizes.

The grid emitted is every sorted triple 2 <= n0 <= n1 <= n2 <= 4 (10 cells),
matching the matrix_q02_n*** rows tracked in ../our_results/cert_summary.csv.  Nothing here
is invented: the bounds are the published Hopcroft-Kerr 1971 and Blaser
1999/2003 formulae carried over verbatim from the matrix-multiplication-search
repo.
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "prev_results.csv")


def hk71(n0, n1, n2):
    if n0 == 2 and n1 == 2 and n2 >= 2:
        return ((7 * n2 + 1) // 2, -1971, "hk71")
    if n0 == 3 and n1 == 2 and n2 == 3:
        return (15, -1971, "hk71")
    if n0 == 3 and n1 == 2 and n2 == 4:
        return (19 - 1e-3, -1971, "hk71")  # no proof
    return (0, 0, "")


def blaser1999lower(l, m, n):
    if n >= l and l >= 2:
        return (l * m + m * n + l - m + n - 3, -1999, "blaser1999lower")
    return (0, 0, "")


def blaser1999fivehalf(n0, n1, n2):
    if n0 == n1 == n2 and n0 >= 3:
        return ((5 * n0**2 + 1) // 2 - 3 * n0, -1999, "blaser1999fivehalf")
    return (0, 0, "")


def blaser2003complexity(n, m, n2):
    if n == n2 and m >= n and n >= 3:
        return (2 * m * n + 2 * n - m - 2, -2003, "blaser2003complexity")
    return (0, 0, "")


def rank_lower_bound_n012(n0, n1, n2):
    return max(
        hk71(n0, n1, n2),
        blaser1999lower(n0, n1, n2),
        blaser1999fivehalf(n0, n1, n2),
        blaser2003complexity(n0, n1, n2),
    )


def rank_lower_bound(n0, n1, n2):
    return max(
        rank_lower_bound_n012(n0, n1, n2),
        rank_lower_bound_n012(n0, n2, n1),
        rank_lower_bound_n012(n1, n0, n2),
        rank_lower_bound_n012(n1, n2, n0),
        rank_lower_bound_n012(n2, n0, n1),
        rank_lower_bound_n012(n2, n1, n0),
    )


def main():
    rows = []
    for n0 in range(2, 5):
        for n1 in range(n0, 5):
            for n2 in range(n1, 5):
                value, _negyear, ref = rank_lower_bound(n0, n1, n2)
                rows.append(
                    {
                        "problem": f"matrix_q02_n{n0}{n1}{n2}",
                        "n0": n0,
                        "n1": n1,
                        "n2": n2,
                        "lb": value,
                        "lb_ref": ref,
                    }
                )
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["problem", "n0", "n1", "n2", "lb", "lb_ref"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
