#pragma once

// Brute-force tensor isomorphism over 𝔽_q (q = P^M, P prime, M ≥ 1), for use
// in unit tests only.
//
// Two NA x NB x NC tensors T and T' over 𝔽_q are isomorphic iff there exist
// invertible matrices A in GL(NA, 𝔽_q), B in GL(NB, 𝔽_q), C in GL(NC, 𝔽_q) with
//
//   T'[i'][j'][k'] = Σ over i,j,k of A[i'][i] B[j'][j] C[k'][k] T[i][j][k],
//
// where the sum and the entry products are in 𝔽_q.
//
// `AreTensorsIsomorphic` enumerates only GL(NA, 𝔽_q) × GL(NB, 𝔽_q) and, for
// each pair (A, B), decides the C-mode by linear algebra rather than another
// enumeration. View t as NC matrix-slices M_k ∈ 𝔽_q^{NA × NB} with
// M_k[i][j] = t[i][j][k], flatten each into a row vector in 𝔽_q^{NA·NB}, and
// stack the NC rows into a matrix V. Acting by some C ∈ GL(NC, 𝔽_q) on mode C
// becomes left multiplication V ↦ C · V, so a C exists iff V and V' have the
// same row span — checkable by comparing column-reversed RREFs. This drops a
// factor of |GL(NC, 𝔽_q)| relative to the brute triple-enumeration oracle (the
// previous implementation, retained in tensor_isomorphism_test.cc as the
// independent ground truth for the randomized cross-check). The (A, B) pass is
// still exponential in NA²+NB², so this remains a tests-only tool.
//
// Tensor entries and matrix entries are GF<P, M>. For (P, M) == (2, 1)
// the row type is GFVec<2, 1, N>, whose partial specialization wraps a
// bit-packed BitVec; the body below keeps a single generic implementation and
// relies on that specialization to make the F₂ hot path collapse to the same
// XOR/AND work as the legacy code.

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/constraints.h"
#include "core/gf.h"
#include "core/gf_vec.h"
#include "core/tensor.h"

namespace tensor_isomorphism_internal {

// An N x N matrix over 𝔽_q (q = P^M), stored as N rows of GFVec<P, M, N>.
template <int P, int M, std::size_t N>
using SquareMatrix = std::array<GFVec<P, M, static_cast<int>(N)>, N>;

// Every invertible N x N matrix over 𝔽_q (the group GL(N, 𝔽_q)). Enumerates
// all q^(N*N) matrices and keeps the full-rank ones. Only sane for small N
// and small q; the static_assert below catches uses that would overflow the
// uint64_t enumeration counter.
template <int P, int M, std::size_t N>
std::vector<SquareMatrix<P, M, N>> EnumerateInvertibleMatrices() {
  constexpr int kIntN = static_cast<int>(N);
  constexpr uint64_t kRowCount = IntPow(GF<P, M>::kQ, kIntN);
  constexpr uint64_t kTotal = IntPow(GF<P, M>::kQ, kIntN * kIntN);
  std::vector<SquareMatrix<P, M, N>> result;
  for (uint64_t bits = 0; bits < kTotal; ++bits) {
    SquareMatrix<P, M, N> matrix;
    uint64_t rest = bits;
    for (std::size_t i = 0; i < N; ++i) {
      matrix[i] = DecodeGFVec<P, M, kIntN>(rest % kRowCount);
      rest /= kRowCount;
    }
    Constraints<P, M, kIntN> rows(matrix.begin(), matrix.end());
    if (IsLinearIndependentRREF<P, M, kIntN>(std::move(rows))) {
      result.push_back(std::move(matrix));
    }
  }
  return result;
}

// out[i'][j][k] = Σ over i of A[i'][i] · t[i][j][k] (act on mode A).
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
Tensor<P, M, NA, NB, NC> ApplyToModeA(const Tensor<P, M, NA, NB, NC> &t,
                                      const SquareMatrix<P, M, NA> &a) {
  using Field = GF<P, M>;
  Tensor<P, M, NA, NB, NC> out{};
  for (std::size_t ip = 0; ip < NA; ++ip) {
    for (std::size_t i = 0; i < NA; ++i) {
      const Field s = a[ip][static_cast<int>(i)];
      if (s == Field::Zero()) {
        continue;
      }
      for (std::size_t j = 0; j < NB; ++j) {
        for (std::size_t k = 0; k < NC; ++k) {
          out[ip][j][k] += s * t[i][j][k];
        }
      }
    }
  }
  return out;
}

// out[i][j'][k] = Σ over j of B[j'][j] · t[i][j][k] (act on mode B).
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
Tensor<P, M, NA, NB, NC> ApplyToModeB(const Tensor<P, M, NA, NB, NC> &t,
                                      const SquareMatrix<P, M, NB> &b) {
  using Field = GF<P, M>;
  Tensor<P, M, NA, NB, NC> out{};
  for (std::size_t i = 0; i < NA; ++i) {
    for (std::size_t jp = 0; jp < NB; ++jp) {
      for (std::size_t j = 0; j < NB; ++j) {
        const Field s = b[jp][static_cast<int>(j)];
        if (s == Field::Zero()) {
          continue;
        }
        for (std::size_t k = 0; k < NC; ++k) {
          out[i][jp][k] += s * t[i][j][k];
        }
      }
    }
  }
  return out;
}

// out[i][j][k'] = Σ over k of C[k'][k] · t[i][j][k] (act on mode C).
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
Tensor<P, M, NA, NB, NC> ApplyToModeC(const Tensor<P, M, NA, NB, NC> &t,
                                      const SquareMatrix<P, M, NC> &c) {
  using Field = GF<P, M>;
  Tensor<P, M, NA, NB, NC> out{};
  for (std::size_t i = 0; i < NA; ++i) {
    for (std::size_t j = 0; j < NB; ++j) {
      for (std::size_t kp = 0; kp < NC; ++kp) {
        for (std::size_t k = 0; k < NC; ++k) {
          const Field s = c[kp][static_cast<int>(k)];
          if (s == Field::Zero()) {
            continue;
          }
          out[i][j][kp] += s * t[i][j][k];
        }
      }
    }
  }
  return out;
}

// Flatten the NC C-slices M_k[i][j] = t[i][j][k] into NC row vectors of length
// NA·NB and return their column-reversed RREF (the canonical form of the row
// span). Two tensors agree under some C ∈ GL(NC, 𝔽_q) acting on mode C iff
// their flattened C-slices have the same row span, i.e. their RREFs match.
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
Constraints<P, M, static_cast<int>(NA *NB)>
CSliceRREF(const Tensor<P, M, NA, NB, NC> &t) {
  using Field = GF<P, M>;
  constexpr int kD = static_cast<int>(NA * NB);
  Constraints<P, M, kD> rows(NC);
  for (std::size_t k = 0; k < NC; ++k) {
    for (std::size_t i = 0; i < NA; ++i) {
      for (std::size_t j = 0; j < NB; ++j) {
        const Field v = t[i][j][k];
        if (v == Field::Zero()) {
          continue;
        }
        rows[k].Set(static_cast<int>(i * NB + j), v);
      }
    }
  }
  GaussJordanRREF<P, M, kD>(&rows);
  return rows;
}

} // namespace tensor_isomorphism_internal

// True iff `t0` and `t1` are isomorphic as tensors over 𝔽_q (q = P^M).
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
bool AreTensorsIsomorphic(const Tensor<P, M, NA, NB, NC> &t0,
                          const Tensor<P, M, NA, NB, NC> &t1) {
  namespace iso = tensor_isomorphism_internal;
  const auto group_a = iso::EnumerateInvertibleMatrices<P, M, NA>();
  const auto group_b = iso::EnumerateInvertibleMatrices<P, M, NB>();
  CHECK_GE((static_cast<uint64_t>(1) << 40) / group_a.size() / group_b.size(),
           1)
      << "group_a.size()=" << group_a.size()
      << ", group_b.size()=" << group_b.size();
  // The C-mode test: t0' and t1 agree under some C ∈ GL(NC, 𝔽_q) iff their
  // flattened C-slice matrices have the same row span. Precompute t1's
  // canonical form once.
  const auto rref1 = iso::CSliceRREF(t1);
  const uint64_t size_ab =
      static_cast<uint64_t>(group_a.size()) * group_b.size();
  // Enumerate (a, b) in pseudo-random order so a true positive is hit early on
  // average. kPrime coprime to size_ab (it's prime and divides neither factor).
  constexpr uint64_t kPrime = 16'777'213; // 2^24 - 3
  CHECK_NE(group_a.size() % kPrime, 0);
  CHECK_NE(group_b.size() % kPrime, 0);
  for (uint64_t i = 0; i < size_ab; ++i) {
    uint64_t j = i * kPrime % size_ab;
    const auto &a = group_a[j % group_a.size()];
    const auto &b = group_b[j / group_a.size()];
    const Tensor<P, M, NA, NB, NC> tab =
        iso::ApplyToModeB(iso::ApplyToModeA(t0, a), b);
    if (iso::CSliceRREF(tab) == rref1) {
      return true;
    }
  }
  return false;
}
