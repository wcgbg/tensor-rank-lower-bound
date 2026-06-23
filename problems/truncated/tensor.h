#pragma once

#include <cstddef>

#include "core/tensor.h"

namespace truncated {

template <int P, int M, std::size_t N>
Tensor<P, M, N, N, N> BuildMulTensor() {
  // The tensor of the truncated multiplication of two (N-1)-degree polynomials
  // in F_q[x]/x^N: c_k = sum_{i+j == k} a_i*b_j keeping only degrees 0..N-1, so
  // T[i][j][k] = 1 iff i+j == k < N.
  Tensor<P, M, N, N, N> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N - i; ++j) {
      tensor[i][j][i + j] = 1;
    }
  }
  return tensor;
}

} // namespace truncated
