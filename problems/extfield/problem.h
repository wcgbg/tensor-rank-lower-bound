#pragma once

// Extension-field multiplication 𝔽_{q^n} over 𝔽_q with q = P^M.
// i.e. multiplication in 𝔽_q[x] / p(x), p the lex-smallest irreducible of
// degree N over 𝔽_q (see irreducibles.h).
//
//   T_ext(N) : 𝔽_{q^n} × 𝔽_{q^n} → 𝔽_{q^n}
//   (f, g) ↦ f · g mod p(x)
//
// As a tensor:
//   T_ext = Σ_{i, j < N} a_i ⊗ b_j ⊗ c_{r(i,j)},  x^{i+j} = Σ_k r(i,j)_k x^k.

#include <format>
#include <string>

#include "core/bit_vec.h"
#include "core/constraints.h"
#include "core/gf_vec.h" // for IntPow
#include "core/tensor.h"
#include "problems/extfield/irreducibles.h"
#include "problems/extfield/symmetry.h"
#include "problems/extfield/tensor.h"

namespace extfield {

template <int P, int M, int N> struct Problem {
  static_assert(N >= 1);
  static_assert(N <= 16); // backtracking ceiling for the M = 1 path
  static_assert(M >= 1);

  static constexpr int kNA = N;
  static constexpr int kNB = N;
  static constexpr int kNC = N;
  static constexpr int kP = P;
  static constexpr int kM = M;
  static constexpr int kQ = static_cast<int>(IntPow(P, M));

  using Vec = GFVec<P, M, N>;
  using SymmetryGroup = extfield::SymmetryGroup<P, M, N>;

  static Tensor<P, M, kNA, kNB, kNC> MakeTensor() {
    if constexpr (P == 2 && M == 1) {
      return BuildMulTensor<N>(IrreduciblePolyBits<N>());
    } else {
      // General base field 𝔽_q = 𝔽_{P^M} (subsumes P > 2, M = 1, where
      // GF<P, 1> is the same mod-P arithmetic as the legacy uint8_t path).
      return BuildMulTensor<P, M, N>(IrreduciblePolyCoeffs<P, M, N>());
    }
  }

  static std::string Name() {
    return std::format("extfield_q{:02}_n{:02}", kQ, N);
  }
};

} // namespace extfield
