#pragma once

// Truncated (short) product, i.e. multiplication in 𝔽_q[x]/(xᴺ) with q = P^M:
// the low N coefficients of the product of two (N−1)-degree polynomials.
//
//   T_trunc(N) : 𝔽_q[x]/(xᴺ) × 𝔽_q[x]/(xᴺ) → 𝔽_q[x]/(xᴺ)
//   (f, g) ↦ f · g mod xᴺ
//
// As a tensor:
//   T_trunc = Σ_{i+j < N} a_i ⊗ b_j ⊗ c_{i+j}.

#include <format>
#include <string>

#include "core/constraints.h"
#include "core/gf_vec.h" // for IntPow
#include "core/tensor.h"
#include "problems/truncated/symmetry.h"
#include "problems/truncated/tensor.h"

namespace truncated {

template <int P, int M, int N> struct Problem {
  static_assert(N >= 1);
  static_assert(N <= 32); // backtracking ceiling for the M = 1 path
  static_assert(M >= 1);

  static constexpr int kNA = N;
  static constexpr int kNB = N;
  static constexpr int kNC = N;
  static constexpr int kP = P;
  static constexpr int kM = M;
  static constexpr int kQ = static_cast<int>(IntPow(P, M));

  using Vec = GFVec<P, M, N>;
  using SymmetryGroup = truncated::SymmetryGroup<P, M, N>;

  static Tensor<P, M, kNA, kNB, kNC> MakeTensor() {
    return BuildMulTensor<P, M, N>();
  }

  static std::string Name() {
    return std::format("truncated_q{:02}_n{:02}", kQ, N);
  }
};

} // namespace truncated
