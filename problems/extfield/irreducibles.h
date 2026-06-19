#pragma once

// Lex-smallest irreducible polynomial of degree N over 𝔽_P, used to fix the
// 𝔽_{P^N} = 𝔽_P[x]/p(x) representation for extfield.
//
// Two interfaces:
//   - IrreduciblePolyBits<N>(): legacy F_2-only, returns std::bitset<N> of
//     low coefficients (x^0..x^{N-1}); the x^N term is implicit. Consumed by
//     the F_2 fast paths in BuildMulTensor and SymmetryGroup.
//   - IrreduciblePolyCoeffs<P, N>(): general, returns std::array<uint8_t, N>
//     of low coefficients in [0, P). For P=2 it agrees with the bitset above
//     (one byte per bit); for P>2 it searches for the lex-smallest irreducible
//     polynomial at runtime.
//
// "Lex-smallest" = the monic irreducible whose integer encoding
// Σ a_i · P^i (low coefs first) is smallest. The values for P=2 below were
// precomputed by trial division and match standard tables.

#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>

#include "core/gf.h"
#include "core/gf_vec.h"

namespace extfield {

template <std::size_t N> std::bitset<N> IrreduciblePolyBits() {
  static_assert(N >= 1 && N <= 16,
                "no irreducible polynomial tabulated for this degree");
  if constexpr (N == 1)
    return std::bitset<N>(0b0); // x
  else if constexpr (N == 2)
    return std::bitset<N>(0b11); // x^2+x+1
  else if constexpr (N == 3)
    return std::bitset<N>(0b011); // x^3+x+1
  else if constexpr (N == 4)
    return std::bitset<N>(0b0011); // x^4+x+1
  else if constexpr (N == 5)
    return std::bitset<N>(0b00101); // x^5+x^2+1
  else if constexpr (N == 6)
    return std::bitset<N>(0b000011); // x^6+x+1
  else if constexpr (N == 7)
    return std::bitset<N>(0b0000011); // x^7+x+1
  else if constexpr (N == 8)
    return std::bitset<N>(0b00011011); // x^8+x^4+x^3+x+1
  else if constexpr (N == 9)
    return std::bitset<N>(0b000000011); // x^9+x+1
  else if constexpr (N == 10)
    return std::bitset<N>(0b0000001001); // x^10+x^3+1
  else if constexpr (N == 11)
    return std::bitset<N>(0b00000000101); // x^11+x^2+1
  else if constexpr (N == 12)
    return std::bitset<N>(0b000000001001); // x^12+x^3+1
  else if constexpr (N == 13)
    return std::bitset<N>(0b0000000011011); // x^13+x^4+x^3+x+1
  else if constexpr (N == 14)
    return std::bitset<N>(0b00000000100001); // x^14+x^5+1
  else if constexpr (N == 15)
    return std::bitset<N>(0b000000000000011); // x^15+x+1
  else if constexpr (N == 16)
    return std::bitset<N>(0b0000000000101011); // x^16+x^5+x^3+x+1
}

// Polynomial coefficient utility used by the runtime irreducibility check.
// `poly` is a vector of low coefficients (poly[0] = constant term); the
// polynomial is monic if poly.back() == 1.
namespace irreducible_internal {

template <int P> inline uint8_t AddP(uint8_t a, uint8_t b) {
  if constexpr (P == 2)
    return static_cast<uint8_t>(a ^ b);
  const int s = static_cast<int>(a) + static_cast<int>(b);
  return static_cast<uint8_t>(s >= P ? s - P : s);
}

template <int P> inline uint8_t SubP(uint8_t a, uint8_t b) {
  if constexpr (P == 2)
    return static_cast<uint8_t>(a ^ b);
  const int s = static_cast<int>(a) - static_cast<int>(b);
  return static_cast<uint8_t>(s < 0 ? s + P : s);
}

template <int P> inline uint8_t MulP(uint8_t a, uint8_t b) {
  if constexpr (P == 2)
    return static_cast<uint8_t>(a & b);
  return static_cast<uint8_t>((static_cast<int>(a) * static_cast<int>(b)) % P);
}

// Polynomial remainder mod a monic divisor d of degree m, reducing p of degree
// ≤ deg(p). Both poly arrays are low-coef-first; returns the remainder coefs
// (degree < m).
template <int P, std::size_t N>
std::array<uint8_t, N> Mod(std::array<uint8_t, N + 1> p,
                           const std::array<uint8_t, N + 1> &d, int m) {
  for (int i = static_cast<int>(N); i >= m; --i) {
    if (p[i] == 0)
      continue;
    const uint8_t c = p[i]; // d[m] == 1 (monic), so quotient digit is c
    for (int j = 0; j <= m; ++j) {
      p[i - m + j] = SubP<P>(p[i - m + j], MulP<P>(c, d[j]));
    }
  }
  std::array<uint8_t, N> rem{};
  for (int i = 0; i < m; ++i) {
    rem[i] = p[i];
  }
  return rem;
}

// True iff every coefficient of `poly` (length N) is zero.
template <std::size_t N> bool IsZero(const std::array<uint8_t, N> &poly) {
  for (std::size_t i = 0; i < N; ++i) {
    if (poly[i] != 0)
      return false;
  }
  return true;
}

// Convert an integer in [0, P^M) to its base-P digit decomposition (low-digit
// first), padded to length M.
template <int P, std::size_t M>
std::array<uint8_t, M> IntToDigits(unsigned long long i) {
  std::array<uint8_t, M> r{};
  for (std::size_t k = 0; k < M; ++k) {
    r[k] = static_cast<uint8_t>(i % P);
    i /= P;
  }
  return r;
}

// --- 𝔽_q (q = P^M) polynomial helpers for the M ≥ 2 irreducible search ---

// Polynomial remainder mod a monic divisor d of degree m over 𝔽_q, reducing p
// of degree ≤ N. Both poly arrays are low-coef-first (coefficients in
// GF<P, M>); d[m] is implicitly One(). Returns the degree-<m remainder.
template <int P, int M, std::size_t N>
std::array<GF<P, M>, N> ModFq(std::array<GF<P, M>, N + 1> p,
                              const std::array<GF<P, M>, N + 1> &d, int m) {
  using F = GF<P, M>;
  for (int i = static_cast<int>(N); i >= m; --i) {
    if (p[i].value == 0)
      continue;
    const F c = p[i]; // d[m] == One() (monic), so the quotient digit is c
    for (int j = 0; j <= m; ++j) {
      p[i - m + j] = F::Sub(p[i - m + j], F::Mul(c, d[j]));
    }
  }
  std::array<F, N> rem{};
  for (int i = 0; i < m; ++i) {
    rem[i] = p[i];
  }
  return rem;
}

// Lex-smallest monic irreducible of degree N over 𝔽_q = 𝔽_{P^M}, returned as
// its low coefficients in GF<P, M> (the x^N term is implicitly One()).
// Lex order is by Σ cᵢ·qⁱ — i.e. the integer encoding DecodeGFVec uses,
// so iterating an integer counter visits candidates in lex order. Trial
// division by every monic divisor of degree 1..⌊N/2⌋ over 𝔽_q. Runtime (q^N is
// modest for the reachable q ≤ 16).
template <int P, int M, int N>
std::array<GF<P, M>, N> IrreduciblePolyCoeffsFq() {
  using F = GF<P, M>;
  constexpr int kQ = F::kQ;
  const uint64_t limit = IntPow(kQ, N);
  for (uint64_t i = 0; i < limit; ++i) {
    const GFVec<P, M, N> cv = DecodeGFVec<P, M, N>(i);
    std::array<F, N + 1> p{};
    for (int k = 0; k < N; ++k)
      p[k] = cv[k];
    p[N] = F::One();

    bool reducible = false;
    for (int m = 1; m * 2 <= N && !reducible; ++m) {
      const uint64_t divisor_limit = IntPow(kQ, m);
      for (uint64_t d_low = 0; d_low < divisor_limit; ++d_low) {
        const GFVec<P, M, N> dv = DecodeGFVec<P, M, N>(d_low);
        std::array<F, N + 1> d{};
        for (int k = 0; k < m; ++k)
          d[k] = dv[k];
        d[m] = F::One(); // monic
        const auto rem = ModFq<P, M, N>(p, d, m);
        bool is_zero = true;
        for (int t = 0; t < m; ++t) {
          if (rem[t].value != 0) {
            is_zero = false;
            break;
          }
        }
        if (is_zero) {
          reducible = true;
          break;
        }
      }
    }
    if (!reducible) {
      std::array<F, N> r{};
      for (int k = 0; k < N; ++k)
        r[k] = cv[k];
      return r;
    }
  }
  // Unreachable: irreducibles of every degree exist over any finite field.
  return std::array<F, N>{};
}

} // namespace irreducible_internal

// Find the lex-smallest monic irreducible polynomial of degree N over 𝔽_P.
// Returns std::array<uint8_t, N> of low coefficients (x^0..x^{N-1}); the x^N
// term is implicit 1.
template <int P, int N> std::array<uint8_t, N> IrreduciblePolyCoeffs() {
  static_assert(N >= 1);
  static_assert(P >= 2 && P < 256);

  // P=2 path: use the hardcoded table for correctness and zero-cost
  // determinism.
  if constexpr (P == 2) {
    auto bits = IrreduciblePolyBits<N>();
    std::array<uint8_t, N> r{};
    for (std::size_t i = 0; i < N; ++i) {
      r[i] = bits[i] ? 1 : 0;
    }
    return r;
  } else {
    namespace ir = irreducible_internal;
    // Search candidates in lex order: integer i = Σ a_j P^j for j < N.
    const unsigned long long limit = IntPow(P, N);
    for (unsigned long long i = 0; i < limit; ++i) {
      auto coefs = ir::IntToDigits<P, N>(i);
      // Build monic polynomial p of degree N: p[0..N-1] = coefs, p[N] = 1.
      std::array<uint8_t, N + 1> p{};
      for (std::size_t k = 0; k < N; ++k)
        p[k] = coefs[k];
      p[N] = 1;

      // Check irreducibility by trial division by every monic poly of degree
      // 1..floor(N/2).
      bool reducible = false;
      for (int m = 1; m * 2 <= static_cast<int>(N) && !reducible; ++m) {
        const unsigned long long divisor_limit = IntPow(P, m);
        for (unsigned long long d_low = 0; d_low < divisor_limit; ++d_low) {
          std::array<uint8_t, N + 1> d{};
          unsigned long long t = d_low;
          for (int k = 0; k < m; ++k) {
            d[k] = static_cast<uint8_t>(t % P);
            t /= P;
          }
          d[m] = 1; // monic
          const auto rem = ir::Mod<P, N>(p, d, m);
          if (ir::IsZero(rem)) {
            reducible = true;
            break;
          }
        }
      }
      if (!reducible) {
        return coefs;
      }
    }
    // Should never happen for prime P: irreducible polys of every degree
    // exist over a finite field.
    return std::array<uint8_t, N>{};
  }
}

// General (P, M) overload: the lex-smallest monic irreducible of degree N over
// the base field 𝔽_q = 𝔽_{P^M}, returned as low coefficients in GF<P, M>
// (the x^N term is implicit One()). For (P, M) == (2, 1) it reuses the
// hardcoded bitset table (zero-cost, deterministic); otherwise it runs the
// 𝔽_q trial-division search. This is the M ≥ 2 entry point used by
// extfield::Problem::MakeTensor and SymmetryGroup.
template <int P, int M, int N> std::array<GF<P, M>, N> IrreduciblePolyCoeffs() {
  static_assert(N >= 1);
  static_assert(M >= 1);
  using F = GF<P, M>;
  if constexpr (P == 2 && M == 1) {
    auto bits = IrreduciblePolyBits<N>();
    std::array<F, N> r{};
    for (std::size_t i = 0; i < N; ++i) {
      r[i] = F{static_cast<uint8_t>(bits[i] ? 1 : 0)};
    }
    return r;
  } else {
    return irreducible_internal::IrreduciblePolyCoeffsFq<P, M, N>();
  }
}

} // namespace extfield
