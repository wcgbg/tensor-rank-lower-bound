#pragma once

#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>

#include "core/gf.h"
#include "core/tensor.h"

namespace extfield {

// Tensor of multiplication in 𝔽_P[x]/(p(x)) where p is monic of degree N and
// irreducible. Two flavors:
//   - BuildMulTensor<N>(std::bitset<N>): F_2-only, legacy interface,
//     bitset holds low coefficients of p. Returns Tensor<2, 1, N, N, N>.
//   - BuildMulTensor<P, N>(std::array<uint8_t, N>): general F_P (M=1),
//     array holds low coefficients of p (each digit in [0, P)). Returns
//     Tensor<P, 1, N, N, N>.

template <std::size_t N>
Tensor<2, 1, N, N, N> BuildMulTensor(const std::bitset<N> &p) {
  // Precompute red[m] = x^m mod p as a bitset for m in [0, 2N-2].
  std::array<std::bitset<N>, 2 * N - 1> red{};
  for (std::size_t m = 0; m < N; ++m) {
    red[m].set(m);
  }
  for (std::size_t m = N; m <= 2 * N - 2; ++m) {
    const bool overflow = red[m - 1][N - 1];
    red[m] = red[m - 1] << 1;
    if (overflow) {
      red[m] ^= p;
    }
  }
  Tensor<2, 1, N, N, N> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      const std::bitset<N> &r = red[i + j];
      for (std::size_t k = 0; k < N; ++k) {
        tensor[i][j][k] = r[k] ? 1 : 0;
      }
    }
  }
  return tensor;
}

// General-P version. `p_low` holds coefficients (x^0..x^{N-1}) of the
// irreducible polynomial; the x^N term is implicit 1. Tensor entries are in
// [0, P).
template <int P, std::size_t N>
Tensor<P, 1, N, N, N> BuildMulTensor(const std::array<uint8_t, N> &p_low) {
  using Field = GF<P, 1>;
  // Precompute red[m] = x^m mod p as an array of low coefficients for m in
  // [0, 2N-2]. red[N] = x^N mod p = -p_low (per-digit Neg).
  std::array<std::array<Field, N>, 2 * N - 1> red{};
  for (std::size_t m = 0; m < N; ++m) {
    red[m][m] = Field{static_cast<uint8_t>(1)};
  }
  for (std::size_t m = N; m <= 2 * N - 2; ++m) {
    // red[m] = x · red[m-1] mod p: left-shift coefficients; if the top
    // (degree N-1) coefficient was c, subtract c · p_low from the shifted
    // polynomial (i.e., the implicit x^N becomes -p_low).
    const Field top = red[m - 1][N - 1];
    for (std::size_t k = N; k > 0; --k) {
      red[m][k - 1] = (k >= 2) ? red[m - 1][k - 2] : Field::Zero();
    }
    // After the left shift, the would-be x^N coefficient is `top`; substitute
    // x^N → -p_low.
    if (top != Field::Zero()) {
      for (std::size_t k = 0; k < N; ++k) {
        red[m][k] = red[m][k] - top * Field{p_low[k]};
      }
    }
  }
  Tensor<P, 1, N, N, N> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      const auto &r = red[i + j];
      for (std::size_t k = 0; k < N; ++k) {
        tensor[i][j][k] = r[k];
      }
    }
  }
  return tensor;
}

// General (P, M) version over the base field 𝔽_q = 𝔽_{P^M}. `p_low` holds the
// coefficients (x^0..x^{N-1}) of the monic irreducible p of degree N over 𝔽_q;
// the x^N term is implicit One(). Tensor entries are GF<P, M>. This is
// the same x^m mod p reduction as the 𝔽_P body, with the coefficient field
// lifted from GF<P, 1> to GF<P, M> (M ≥ 1; subsumes the array<uint8_t> overload
// above when M == 1, but that one is kept for the existing call sites).
template <int P, int M, std::size_t N>
Tensor<P, M, N, N, N> BuildMulTensor(const std::array<GF<P, M>, N> &p_low) {
  using Field = GF<P, M>;
  // red[m] = x^m mod p as low coefficients, for m in [0, 2N-2].
  std::array<std::array<Field, N>, 2 * N - 1> red{};
  for (std::size_t m = 0; m < N; ++m) {
    red[m][m] = Field::One();
  }
  for (std::size_t m = N; m <= 2 * N - 2; ++m) {
    // red[m] = x · red[m-1] mod p: left-shift coefficients; the implicit x^N
    // becomes -p_low, so subtract (top coefficient) · p_low.
    const Field top = red[m - 1][N - 1];
    for (std::size_t k = N; k > 0; --k) {
      red[m][k - 1] = (k >= 2) ? red[m - 1][k - 2] : Field::Zero();
    }
    if (!(top == Field::Zero())) {
      for (std::size_t k = 0; k < N; ++k) {
        red[m][k] = Field::Sub(red[m][k], Field::Mul(top, p_low[k]));
      }
    }
  }
  Tensor<P, M, N, N, N> tensor{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      const auto &r = red[i + j];
      for (std::size_t k = 0; k < N; ++k) {
        tensor[i][j][k] = r[k];
      }
    }
  }
  return tensor;
}

} // namespace extfield
