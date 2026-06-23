#pragma once

// A-side symmetry group for the truncated (short) product over 𝔽_q with
// q = P^M: multiplication in the local algebra R = 𝔽_q[x]/(xᴺ).
//
//   T_trunc(N) : R × R → R,  T_trunc = Σ_{i+j < N} a_i ⊗ b_j ⊗ c_{i+j}.
//
// truncated is the γ = 0 case of the shared polynomial-quotient symmetry engine
// (R = 𝔽_q[x]/(xᴺ − γ)); see problems/poly_quotient/symmetry.h. Here x is
// nilpotent (xᴺ = 0): the units are {p : p_0 ≠ 0} (so |Store| = q^{N−1}) and the
// ring autos send x ↦ y with yᴺ = 0 and 1, y, …, y^{N−1} a basis (so
// |Query| = M·(q−1)·q^{N−2}, and 1 for N = 1 where x ≡ 0).

#include "problems/poly_quotient/symmetry.h"

namespace truncated {

// xᴺ ≡ 0, so kGamma = 0.
template <int P, int M, int N>
using SymmetryGroup = poly_quotient::SymmetryGroup<P, M, N, 0>;

} // namespace truncated
