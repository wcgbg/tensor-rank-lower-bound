#pragma once

// Galois field 𝔽_q with q = P^M (P prime, M ≥ 1).
//
// GF<P, M> is BOTH:
//   - a value type: one byte of storage, the F_q index in [0, q).
//   - the field's static-ops home: kQ, kGroundIrreducible, Cayley tables, and
//     static Add/Sub/Mul/Neg/Inverse all live on the type itself.
//
// Element bijection: index k ↔ polynomial Σ a_i y^i where a_i = (k / P^i) mod P
// (digit expansion in base P, low coefficient first). So index 0 → 0,
// 1 → 1, P → y, P² → y², etc.
//
//   - M == 1 (partial specialization): the prime field 𝔽_P. Direct mod-P
//     arithmetic (no tables); P == 2 collapses to XOR/AND. Inverse via Fermat.
//   - M ≥ 2 (primary template): 𝔽_{P^M} ≅ 𝔽_P[y] / g(y), where g(y) is the
//     lex-smallest monic irreducible of degree M over 𝔽_P (constexpr trial
//     division). Arithmetic goes through constexpr Cayley tables. q ≤ 16 so
//     each binary table is ≤ 256 bytes and each unary table is ≤ 16 bytes.
//
// Layout: sizeof(GF<P, M>) == 1, trivially copyable, layout-compatible with
// uint8_t. This lets std::array<GF, N> share its byte representation with
// std::array<uint8_t, N> — load-bearing for the memcpy paths in
// core/constraints.h.

#include <array>
#include <cstdint>
#include <string>
#include <type_traits>

#include <ng-log/logging.h>

namespace intpow_detail {
// Non-constexpr on purpose: reaching this call during constant evaluation makes
// the enclosing expression non-constant, so a compile-time IntPow overflow is a
// required diagnostic (build error) with a clear "non-constexpr function"
// message. The runtime overflow path uses LOG(FATAL) instead. Keep it
// non-constexpr.
inline void Overflow() {}
} // namespace intpow_detail

// base^exp over uint64. Overflow is treated as an error:
//   - constant-evaluated (static_assert, constexpr init, template args): a hard
//     compile error, via the non-constexpr intpow_detail::Overflow().
//   - at runtime: LOG(FATAL). No caller may rely on a saturated return value;
//     a caller that legitimately probes a possibly-huge search-space size must
//     bound the exponent itself before calling (see
//     RankLowerBoundForcedProductA, which caps its iteration count without ever
//     overflowing).
constexpr uint64_t IntPow(int base, int exp) {
  uint64_t r = 1;
  for (int i = 0; i < exp; ++i) {
    if (r > (UINT64_MAX / static_cast<uint64_t>(base))) {
      if (std::is_constant_evaluated()) {
        intpow_detail::Overflow();
      } else {
        LOG(FATAL) << "IntPow overflow: " << base << "^" << exp
                   << " does not fit in uint64";
      }
    }
    r *= static_cast<uint64_t>(base);
  }
  return r;
}

namespace gf_internal {

// Compile-time primality check.
constexpr bool IsPrime(int p) {
  if (p < 2)
    return false;
  if (p == 2)
    return true;
  if (p % 2 == 0)
    return false;
  for (int d = 3; d * d <= p; d += 2) {
    if (p % d == 0)
      return false;
  }
  return true;
}

// Scalar add / sub / mul / neg mod P. Used by GF<P, 1> directly and by the
// per-digit / per-coefficient operations that populate the M ≥ 2 tables.
template <int P> constexpr uint8_t AddModP(uint8_t a, uint8_t b) {
  if constexpr (P == 2) {
    return static_cast<uint8_t>(a ^ b);
  } else {
    const int s = static_cast<int>(a) + static_cast<int>(b);
    return static_cast<uint8_t>(s >= P ? s - P : s);
  }
}

template <int P> constexpr uint8_t SubModP(uint8_t a, uint8_t b) {
  if constexpr (P == 2) {
    return static_cast<uint8_t>(a ^ b);
  } else {
    const int s = static_cast<int>(a) - static_cast<int>(b);
    return static_cast<uint8_t>(s < 0 ? s + P : s);
  }
}

template <int P> constexpr uint8_t MulModP(uint8_t a, uint8_t b) {
  if constexpr (P == 2) {
    return static_cast<uint8_t>(a & b);
  } else {
    return static_cast<uint8_t>((static_cast<int>(a) * static_cast<int>(b)) %
                                P);
  }
}

template <int P> constexpr uint8_t NegModP(uint8_t a) {
  if constexpr (P == 2) {
    return a;
  } else {
    return static_cast<uint8_t>(a == 0 ? 0 : P - a);
  }
}

// True iff the monic polynomial of degree deg whose low coefficients are
// `poly[0..deg-1]` (and whose deg-th coefficient is implicitly 1) is divisible
// by *some* monic polynomial of degree 1 .. deg/2 over 𝔽_P.
template <int P, int MaxDeg>
constexpr bool
IsReducibleByLowerDegree(const std::array<uint8_t, MaxDeg + 1> &poly, int deg) {
  for (int m = 1; m * 2 <= deg; ++m) {
    const uint64_t divisor_limit = IntPow(P, m);
    for (uint64_t d_low = 0; d_low < divisor_limit; ++d_low) {
      std::array<uint8_t, MaxDeg + 1> d{};
      uint64_t t = d_low;
      for (int k = 0; k < m; ++k) {
        d[k] = static_cast<uint8_t>(t % P);
        t /= P;
      }
      d[m] = 1;
      std::array<uint8_t, MaxDeg + 1> r = poly;
      for (int i = deg; i >= m; --i) {
        if (r[i] == 0)
          continue;
        const uint8_t c = r[i];
        for (int j = 0; j <= m; ++j) {
          r[i - m + j] = SubModP<P>(r[i - m + j], MulModP<P>(c, d[j]));
        }
      }
      bool is_zero = true;
      for (int i = 0; i < m; ++i) {
        if (r[i] != 0) {
          is_zero = false;
          break;
        }
      }
      if (is_zero) {
        return true;
      }
    }
  }
  return false;
}

// The lex-smallest monic irreducible polynomial of degree M over 𝔽_P, returned
// as its low coefficients (length M; the M-th coefficient is implicitly 1).
template <int P, int M> constexpr std::array<uint8_t, M> GroundIrreducible() {
  static_assert(M >= 1);
  const uint64_t limit = IntPow(P, M);
  for (uint64_t i = 0; i < limit; ++i) {
    std::array<uint8_t, M + 1> p{};
    uint64_t t = i;
    for (int k = 0; k < M; ++k) {
      p[k] = static_cast<uint8_t>(t % P);
      t /= P;
    }
    p[M] = 1;
    if (!IsReducibleByLowerDegree<P, M>(p, M)) {
      std::array<uint8_t, M> r{};
      for (int k = 0; k < M; ++k)
        r[k] = p[k];
      return r;
    }
  }
  return std::array<uint8_t, M>{}; // unreachable for prime P, any M ≥ 1
}

// Decode an F_q index k into its M base-P digits, low-coefficient first.
template <int P, int M>
constexpr std::array<uint8_t, M> IndexToDigits(uint8_t k) {
  std::array<uint8_t, M> d{};
  uint64_t t = k;
  for (int i = 0; i < M; ++i) {
    d[i] = static_cast<uint8_t>(t % P);
    t /= P;
  }
  return d;
}

// Encode M base-P digits (low-coef first) into an F_q index Σ d[i] · P^i.
template <int P, int M>
constexpr uint8_t DigitsToIndex(const std::array<uint8_t, M> &d) {
  uint64_t r = 0;
  for (int i = M - 1; i >= 0; --i) {
    r = r * static_cast<uint64_t>(P) + static_cast<uint64_t>(d[i]);
  }
  return static_cast<uint8_t>(r);
}

// Schoolbook poly multiply mod g(y), on digit-vector inputs (low-coef first).
// Used at compile time to populate GF<P, M>::kMulTable. Callers must pass g
// once (e.g. kGroundIrreducible) — do not recompute GroundIrreducible inside
// the loop.
template <int P, int M>
constexpr std::array<uint8_t, M> MulDigits(const std::array<uint8_t, M> &a,
                                           const std::array<uint8_t, M> &b,
                                           const std::array<uint8_t, M> &g) {
  std::array<uint8_t, 2 * M> prod{};
  for (int i = 0; i < M; ++i) {
    if (a[i] == 0)
      continue;
    for (int j = 0; j < M; ++j) {
      if (b[j] == 0)
        continue;
      prod[i + j] = AddModP<P>(prod[i + j], MulModP<P>(a[i], b[j]));
    }
  }
  for (int k = 2 * M - 2; k >= M; --k) {
    const uint8_t c = prod[k];
    if (c == 0)
      continue;
    for (int j = 0; j < M; ++j) {
      prod[k - M + j] = SubModP<P>(prod[k - M + j], MulModP<P>(c, g[j]));
    }
    prod[k] = 0;
  }
  std::array<uint8_t, M> r{};
  for (int i = 0; i < M; ++i)
    r[i] = prod[i];
  return r;
}

} // namespace gf_internal

// Primary template: M ≥ 2 (the M = 1 specialization below supersedes this body
// for the prime-field hot path). Value-type + static-ops in one struct.
template <int P, int M> struct GF {
  static_assert(P >= 2 && P < 256, "P must satisfy 2 <= P < 256");
  static_assert(M >= 1, "M must be at least 1");
  static_assert(gf_internal::IsPrime(P), "P must be prime");
  static_assert(IntPow(P, M) <= 16,
                "q = P^M must be <= 16 for the lookup-table GF");

  uint8_t value = 0;
  static constexpr int kQ = static_cast<int>(IntPow(P, M));
  static constexpr std::array<uint8_t, M> kGroundIrreducible =
      gf_internal::GroundIrreducible<P, M>();

  // Ctors / conversion. Implicit `GF(uint8_t)` so `gf = 1`, `gf = 0` are
  // idiomatic; the explicit-to-uint8_t conversion is for the rare callers
  // that need the raw index (e.g., sparse-string emission, hashing).
  constexpr GF() = default;
  constexpr GF(uint8_t v) : value(v) {}
  constexpr explicit operator uint8_t() const { return value; }

  static constexpr GF Zero() { return GF{static_cast<uint8_t>(0)}; }
  static constexpr GF One() { return GF{static_cast<uint8_t>(1)}; }

  constexpr bool operator==(const GF &) const = default;

private:
  using BinaryTable = std::array<std::array<uint8_t, kQ>, kQ>;
  using UnaryTable = std::array<uint8_t, kQ>;

  static constexpr BinaryTable BuildAddTable() {
    BinaryTable t{};
    for (int a = 0; a < kQ; ++a) {
      const auto da = gf_internal::IndexToDigits<P, M>(static_cast<uint8_t>(a));
      for (int b = 0; b < kQ; ++b) {
        const auto db =
            gf_internal::IndexToDigits<P, M>(static_cast<uint8_t>(b));
        std::array<uint8_t, M> ds{};
        for (int i = 0; i < M; ++i) {
          ds[i] = gf_internal::AddModP<P>(da[i], db[i]);
        }
        t[a][b] = gf_internal::DigitsToIndex<P, M>(ds);
      }
    }
    return t;
  }

  static constexpr BinaryTable
  BuildMulTable(const std::array<uint8_t, M> &ground) {
    BinaryTable t{};
    for (int a = 0; a < kQ; ++a) {
      const auto da = gf_internal::IndexToDigits<P, M>(static_cast<uint8_t>(a));
      for (int b = 0; b < kQ; ++b) {
        const auto db =
            gf_internal::IndexToDigits<P, M>(static_cast<uint8_t>(b));
        t[a][b] = gf_internal::DigitsToIndex<P, M>(
            gf_internal::MulDigits<P, M>(da, db, ground));
      }
    }
    return t;
  }

  static constexpr UnaryTable BuildInvTable(const BinaryTable &mul) {
    UnaryTable t{}; // t[0] = 0; Inverse(0) is ill-defined but returns 0.
    for (int a = 1; a < kQ; ++a) {
      for (int b = 1; b < kQ; ++b) {
        if (mul[a][b] == 1) {
          t[a] = static_cast<uint8_t>(b);
          break;
        }
      }
    }
    return t;
  }

public:
  static constexpr BinaryTable kAddTable = BuildAddTable();
  static constexpr BinaryTable kMulTable = BuildMulTable(kGroundIrreducible);
  static constexpr UnaryTable kInvTable = BuildInvTable(kMulTable);

  // Static ops — take and return GF.
  static constexpr GF Add(GF a, GF b) {
    return GF{kAddTable[a.value][b.value]};
  }
  static constexpr GF Mul(GF a, GF b) {
    return GF{kMulTable[a.value][b.value]};
  }
  static constexpr GF Inverse(GF a) { return GF{kInvTable[a.value]}; }

  // a^exp via binary exponentiation. Pow(0, 0) = 1 by convention.
  static constexpr GF Pow(GF a, int exp) {
    GF result = One();
    GF base = a;
    while (exp > 0) {
      if (exp & 1) {
        result = Mul(result, base);
      }
      base = Mul(base, base);
      exp >>= 1;
    }
    return result;
  }

  // Additive inverse. Char-2 fields (F_4, F_8, F_16): identity. F_9: hard-
  // coded permutation. Anything else fires the static_assert below.
  static constexpr GF Neg(GF a) {
    if constexpr (P == 2) {
      return a;
    } else if constexpr (P == 3 && M == 2) {
      constexpr uint8_t kF9Neg[9] = {0, 2, 1, 6, 8, 7, 3, 5, 4};
      return GF{kF9Neg[a.value]};
    } else {
      static_assert(P == 2 || (P == 3 && M == 2),
                    "GF::Neg: only char-2 fields and F_9 are hard-coded; add "
                    "a branch when extending the supported (P, M) set.");
    }
  }

  static constexpr GF Sub(GF a, GF b) { return Add(a, Neg(b)); }

  // Operators (delegate to statics).
  friend constexpr GF operator+(GF a, GF b) { return Add(a, b); }
  friend constexpr GF operator-(GF a, GF b) { return Sub(a, b); }
  friend constexpr GF operator*(GF a, GF b) { return Mul(a, b); }
  friend constexpr GF operator-(GF a) { return Neg(a); }
  constexpr GF Inverse() const { return Inverse(*this); }

  constexpr GF &operator+=(GF b) {
    value = Add(*this, b).value;
    return *this;
  }
  constexpr GF &operator-=(GF b) {
    value = Sub(*this, b).value;
    return *this;
  }
  constexpr GF &operator*=(GF b) {
    value = Mul(*this, b).value;
    return *this;
  }

  // Concatenated base-P digits (low coefficient first), one char per digit.
  // q = P^M ≤ 16 with M ≥ 2 forces P ∈ {2, 3}, so every digit fits in '0'..'2'.
  std::string ToString() const {
    const auto digits = gf_internal::IndexToDigits<P, M>(value);
    std::string s;
    s.reserve(M);
    for (int i = 0; i < M; ++i) {
      s.push_back(static_cast<char>('0' + digits[i]));
    }
    return s;
  }
};

// Partial specialization for M = 1: prime field 𝔽_P. Direct mod-P arithmetic;
// byte-identical codegen to the legacy FieldOps<P, 1>.
template <int P> struct GF<P, 1> {
  static_assert(P >= 2 && P < 256, "P must satisfy 2 <= P < 256");
  static_assert(gf_internal::IsPrime(P), "P must be prime");

  uint8_t value = 0;
  static constexpr int kQ = P;

  constexpr GF() = default;
  constexpr GF(uint8_t v) : value(v) {}
  constexpr explicit operator uint8_t() const { return value; }

  static constexpr GF Zero() { return GF{static_cast<uint8_t>(0)}; }
  static constexpr GF One() { return GF{static_cast<uint8_t>(1)}; }

  constexpr bool operator==(const GF &) const = default;

  static constexpr GF Add(GF a, GF b) {
    return GF{gf_internal::AddModP<P>(a.value, b.value)};
  }
  static constexpr GF Sub(GF a, GF b) {
    return GF{gf_internal::SubModP<P>(a.value, b.value)};
  }
  static constexpr GF Mul(GF a, GF b) {
    return GF{gf_internal::MulModP<P>(a.value, b.value)};
  }
  static constexpr GF Neg(GF a) { return GF{gf_internal::NegModP<P>(a.value)}; }

  // Multiplicative inverse via Fermat: a^(P−2) = a⁻¹. Inverse(0) returns 0
  // (ill-defined but does not produce surprising values; callers must not
  // invert 0).
  static constexpr GF Inverse(GF a) {
    if (a.value == 0) {
      return Zero();
    }
    if constexpr (P == 2) {
      return One();
    } else {
      int result = 1;
      int base = a.value;
      int e = P - 2;
      while (e > 0) {
        if (e & 1) {
          result = (result * base) % P;
        }
        base = (base * base) % P;
        e >>= 1;
      }
      return GF{static_cast<uint8_t>(result)};
    }
  }

  // a^exp via binary exponentiation. Pow(0, 0) = 1 by convention.
  static constexpr GF Pow(GF a, int exp) {
    GF result = One();
    GF base = a;
    while (exp > 0) {
      if (exp & 1) {
        result = Mul(result, base);
      }
      base = Mul(base, base);
      exp >>= 1;
    }
    return result;
  }

  friend constexpr GF operator+(GF a, GF b) { return Add(a, b); }
  friend constexpr GF operator-(GF a, GF b) { return Sub(a, b); }
  friend constexpr GF operator*(GF a, GF b) { return Mul(a, b); }
  friend constexpr GF operator-(GF a) { return Neg(a); }
  constexpr GF Inverse() const { return Inverse(*this); }

  constexpr GF &operator+=(GF b) {
    value = Add(*this, b).value;
    return *this;
  }
  constexpr GF &operator-=(GF b) {
    value = Sub(*this, b).value;
    return *this;
  }
  constexpr GF &operator*=(GF b) {
    value = Mul(*this, b).value;
    return *this;
  }

  // Decimal representation of the F_P index — multi-character for P > 10
  // (e.g., GF<13, 1>{12} → "12") so values stay unambiguous.
  std::string ToString() const {
    return std::to_string(static_cast<int>(value));
  }
};

// Layout guarantees the rest of the codebase relies on.
static_assert(sizeof(GF<2, 1>) == 1);
static_assert(sizeof(GF<3, 1>) == 1);
static_assert(sizeof(GF<2, 2>) == 1);
static_assert(sizeof(GF<3, 2>) == 1);
static_assert(std::is_trivially_copyable_v<GF<2, 1>>);
static_assert(std::is_trivially_copyable_v<GF<3, 2>>);
