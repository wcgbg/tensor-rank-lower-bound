#pragma once

// A-side symmetry group for cyclic polynomial multiplication over 𝔽_q with
// q = P^M: multiplication in R = 𝔽_q[x]/(xᴺ−1).
//
//   T_cyc(N) : R × R → R,  T_cyc = Σ_{i,j < N} a_i ⊗ b_j ⊗ c_{(i+j) mod N}.
//
// cyclic is the γ = 1 case of the shared polynomial-quotient symmetry engine
// (R = 𝔽_q[x]/(xᴺ − γ)); see problems/poly_quotient/symmetry.h for the full
// derivation of the group G = (R^* ⋊ Aut(R) ⋊ Gal)/𝔽_q^* and its meet-in-the-
// middle split Store = R^*/𝔽_q^*, Query = Aut(R) ⋊ Gal.

#include "problems/poly_quotient/symmetry.h"

namespace cyclic {

// xᴺ ≡ 1, so kGamma = 1.
template <int P, int M, int N>
using SymmetryGroup = poly_quotient::SymmetryGroup<P, M, N, 1>;

} // namespace cyclic
