#pragma once

#include <cstddef>
#include <cstdint>

#include "core/gf.h"
#include "core/tensor.h"

namespace negacyclic {

template <int P, int M, std::size_t N>
Tensor<P, M, N, N, N> BuildMulTensor() {
  // The tensor of negacyclic multiplication of two (N-1)-degree polynomials in
  // F_q[x]/(x^N + 1): x^N ≡ -1, so c_k = sum_{i+j == k} a_i*b_j
  //                                     - sum_{i+j == k+N} a_i*b_j.
  // Hence T[i][j][i+j] = 1 when i+j < N, and T[i][j][i+j-N] = -1 when i+j >= N.
  // Over characteristic 2, -1 = 1 and this coincides with the cyclic tensor
  // (x^N + 1 = x^N - 1).
  using F = GF<P, M>;
  const F one = F::One();
  const F neg_one = F{static_cast<uint8_t>(P - 1)}; // -1 ∈ F_p ⊆ F_q
  Tensor<P, M, N, N, N> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      const std::size_t s = i + j;
      if (s < N) {
        tensor[i][j][s] = one;
      } else {
        tensor[i][j][s - N] = neg_one;
      }
    }
  }
  return tensor;
}

} // namespace negacyclic
