#pragma once

#include <cstddef>

#include "core/tensor.h"

namespace cyclic {

template <int P, int M, std::size_t N> Tensor<P, M, N, N, N> BuildMulTensor() {
  // The tensor of the cyclic multiplication of two (N-1)-degree polynomials in
  // F_q[x]/(x^N-1): c_k = sum_{i+j == k (mod N)} a_i*b_j, so T[i][j][k] = 1 iff
  // (i+j) mod N == k.
  Tensor<P, M, N, N, N> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      tensor[i][j][(i + j) % N] = 1;
    }
  }
  return tensor;
}

} // namespace cyclic
