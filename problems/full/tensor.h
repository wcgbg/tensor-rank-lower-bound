#pragma once

#include <cstddef>

#include "core/tensor.h"

namespace full {

template <int P, int M, std::size_t N>
Tensor<P, M, N, N, 2 * N - 1> BuildMulTensor() {
  // The tensor of the multiplication of two (N-1)-degree polynomials, i.e. the
  // convolution: c_k = sum_{i+j=k} a_i*b_j, so T[i][j][k] = 1 iff i+j == k.
  Tensor<P, M, N, N, 2 * N - 1> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      tensor[i][j][i + j] = 1;
    }
  }
  return tensor;
}

} // namespace full
