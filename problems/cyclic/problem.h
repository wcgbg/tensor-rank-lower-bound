#pragma once

// Cyclic convolution, i.e. multiplication in 𝔽_q[x]/(xᴺ−1) with q = P^M.
//
//   T_cyc(N) : 𝔽_q[x]/(xᴺ−1) × 𝔽_q[x]/(xᴺ−1) → 𝔽_q[x]/(xᴺ−1)
//   (f, g) ↦ f · g mod (xᴺ − 1)
//
// As a tensor:
//   T_cyc = Σ_{i, j < N} a_i ⊗ b_j ⊗ c_{(i+j) mod N}.

#include <format>
#include <string>

#include "core/constraints.h"
#include "core/gf_vec.h" // for IntPow
#include "core/tensor.h"
#include "problems/cyclic/symmetry.h"
#include "problems/cyclic/tensor.h"

namespace cyclic {

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
  using SymmetryGroup = cyclic::SymmetryGroup<P, M, N>;

  static Tensor<P, M, kNA, kNB, kNC> MakeTensor() {
    return BuildMulTensor<P, M, N>();
  }

  static std::string Name() {
    return std::format("cyclic_q{:02}_n{:02}", kQ, N);
  }
};

} // namespace cyclic
