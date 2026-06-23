#pragma once

// A-side symmetry group for negacyclic polynomial multiplication over 𝔽_q with
// q = P^M: multiplication in R = 𝔽_q[x]/(xᴺ + 1).
//
//   T_neg(N) : R × R → R,
//   T_neg = Σ_{i+j < N} a_i ⊗ b_j ⊗ c_{i+j}
//         − Σ_{i+j ≥ N} a_i ⊗ b_j ⊗ c_{i+j−N}.
//
// negacyclic is the γ = −1 case of the shared polynomial-quotient symmetry
// engine (R = 𝔽_q[x]/(xᴺ − γ)); see problems/poly_quotient/symmetry.h. The
// reduction constant is the 𝔽_p element −1, whose GFVec index is P−1. Over
// characteristic 2, P−1 = 1, so this is *identical* to cyclic (xᴺ + 1 = xᴺ − 1).
// The unit/auto enumeration follows the factorisation of xᴺ + 1 (the q-
// cyclotomic cosets of the odd residues mod 2N for p odd).

#include "problems/poly_quotient/symmetry.h"

namespace negacyclic {

// xᴺ ≡ −1, so kGamma = P−1 (the 𝔽_p index of −1).
template <int P, int M, int N>
using SymmetryGroup = poly_quotient::SymmetryGroup<P, M, N, P - 1>;

} // namespace negacyclic
