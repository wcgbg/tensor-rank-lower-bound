#!/usr/bin/env python3
"""Generate the lex-smallest irreducible polynomials over F2 for irreducibles.h.

A polynomial over F2 is encoded as an integer: bit d is the coefficient of x^d.
For each degree N we search candidates x^N + v (v = 0, 1, 2, ...) in increasing
integer value and emit the first irreducible one -- this is the "lex-smallest"
monic irreducible. Only the low bits v (coefficients of x^0..x^(N-1)) go into
irreducibles.h; the x^N term is implicit.

A degree-N polynomial is irreducible over F2 iff no polynomial of degree
1..N/2 divides it (any proper factor contains an irreducible factor of degree
<= N/2). We test that by trial division.

Usage:
    python3 problems/extfield/gen_irreducibles.py [max_degree]   # default 16

The output mirrors the branches in irreducibles.h, so the header can be
regenerated / checked by eye.
"""

import sys


def poly_mod(a: int, b: int) -> int:
    """Remainder of polynomial a divided by polynomial b over F2."""
    db = b.bit_length() - 1
    while a and a.bit_length() - 1 >= db:
        a ^= b << ((a.bit_length() - 1) - db)
    return a


def is_irreducible(p: int, n: int) -> bool:
    """True iff the degree-n polynomial p is irreducible over F2."""
    if n == 1:
        return True  # both x and x+1 are irreducible
    for g in range(2, 1 << (n // 2 + 1)):  # all polys of degree 1..n/2
        if poly_mod(p, g) == 0:
            return False
    return True


def lex_smallest_irreducible(n: int) -> int:
    """Low-bits value v of the lex-smallest monic irreducible of degree n."""
    for v in range(1 << n):
        if is_irreducible((1 << n) | v, n):
            return v
    raise AssertionError(f"no irreducible of degree {n}")  # unreachable


def poly_str(p: int) -> str:
    terms = []
    for d in range(p.bit_length() - 1, -1, -1):
        if (p >> d) & 1:
            terms.append("x^%d" % d if d > 1 else ("x" if d == 1 else "1"))
    return "+".join(terms)


def main() -> None:
    max_degree = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    for n in range(1, max_degree + 1):
        v = lex_smallest_irreducible(n)
        p = (1 << n) | v
        # Matches the form of the branches in irreducibles.h.
        print("  else if constexpr (N == %d)" % n)
        print(
            "    return std::bitset<N>(0b%0*d); // %s"
            % (n, int(bin(v)[2:]), poly_str(p))
        )


if __name__ == "__main__":
    main()
