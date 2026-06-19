#pragma once

// ⟨N0, N1, N2⟩ matrix multiplication over 𝔽_q with q = P^M.
//
//   T_mat(N0,N1,N2) : 𝔽_q^{N0×N1} × 𝔽_q^{N1×N2} → 𝔽_q^{N0×N2}
//   (X, Y) ↦ X · Y
//
// As a tensor (the A/B/C factors flattened row-major; see tensor.h):
//   T_mat = Σ_{i<N0, j<N1, k<N2} x_{ij} ⊗ y_{jk} ⊗ z_{ki},
// so kNA = N0·N1, kNB = N1·N2, kNC = N2·N0.
//
// This is the matrix-multiplication problem family, the third bilinear map
// (besides full and cyclic polynomial multiplication) the framework targets.
// Its A-side symmetry group is the √|G| meet-in-the-middle GL_{N0} × GL_{N1}
// split that core/symmetry.h was designed around — see symmetry.h.

#include <format>
#include <string>

#include "core/constraints.h"
#include "core/gf_vec.h" // for IntPow
#include "core/tensor.h"
#include "problems/matrix/symmetry.h"
#include "problems/matrix/tensor.h"

namespace matrix {

template <int P, int M, int N0, int N1, int N2> struct Problem {
  static_assert(N0 >= 1 && N1 >= 1 && N2 >= 1);
  static_assert(N0 * N1 <= 32); // kNA backtracking ceiling for the M = 1 path
  static_assert(M >= 1);
  // The matrix family is 𝔽₂-only: its symmetry group
  // (problems/matrix/symmetry.h) uses 𝔽₂ lookup tables, and the paper's matrix
  // results and the flip-graph upper bound are all over 𝔽₂.
  static_assert(P == 2 && M == 1,
                "matrix::Problem is implemented for 𝔽₂ (P=2, M=1) only");

  static constexpr int kNA = N0 * N1;
  static constexpr int kNB = N1 * N2;
  static constexpr int kNC = N2 * N0;
  static constexpr int kP = P;
  static constexpr int kM = M;
  static constexpr int kQ = static_cast<int>(IntPow(P, M));

  // The three matrix-multiplication dimensions, exposed for tooling.
  static constexpr int kN0 = N0;
  static constexpr int kN1 = N1;
  static constexpr int kN2 = N2;

  using Vec = GFVec<P, M, kNA>;
  using SymmetryGroup = matrix::SymmetryGroup<P, M, N0, N1, N2>;

  static Tensor<P, M, kNA, kNB, kNC> MakeTensor() {
    return BuildMulTensor<P, M, N0, N1, N2>();
  }

  // <family>_q<QQ>_n<N0><N1><N2>, e.g. matrix_q02_n333. All feasible formats
  // have single-digit dimensions (mirroring the source repo's rmms_n333).
  static std::string Name() {
    return std::format("matrix_q{:02}_n{}{}{}", kQ, N0, N1, N2);
  }
};

} // namespace matrix
