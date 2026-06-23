#include "problems/extfield/tensor.h"

#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>

#include <gtest/gtest.h>

#include "core/gf.h"
#include "core/gf_vec.h"
#include "core/tensor.h"
#include "core/tensor_utils.h"
#include "problems/cyclic/tensor.h"
#include "problems/extfield/irreducibles.h"
#include "problems/truncated/tensor.h"

namespace extfield {

// BuildMulTensor over F2[x]/(x^2+x+1): x^2 = x+1.
TEST(BuildMulTensorTest, GF4) {
  std::bitset<2> p;
  p.set(0); // x^0
  p.set(1); // x^1
  // a1*b1 = x^2 = x+1 -> c0 + c1.
  EXPECT_EQ((TensorToSparseString<2, 1, 2, 2, 2>(BuildMulTensor<2>(p))),
            "a0*b0*c0 + a0*b1*c1 + a1*b0*c1 + a1*b1*c0 + a1*b1*c1");
}

// BuildMulTensor over F2[x]/(x^3+x+1): x^3 = x+1.
TEST(BuildMulTensorTest, GF8) {
  std::bitset<3> p;
  p.set(0); // x^0
  p.set(1); // x^1
  const auto tensor = BuildMulTensor<3>(p);
  // Expected x^(i+j) mod (x^3+x+1):
  //   x^0=1, x^1=x, x^2=x^2, x^3=x+1, x^4=x^2+x.
  Tensor<2, 1, 3, 3, 3> expected{};
  expected[0][0][0] = 1;                     // 1
  expected[0][1][1] = 1;                     // x
  expected[0][2][2] = 1;                     // x^2
  expected[1][0][1] = 1;                     // x
  expected[1][1][2] = 1;                     // x^2
  expected[1][2][0] = expected[1][2][1] = 1; // x^3 = x+1
  expected[2][0][2] = 1;                     // x^2
  expected[2][1][0] = expected[2][1][1] = 1; // x^3 = x+1
  expected[2][2][1] = expected[2][2][2] = 1; // x^4 = x^2+x
  EXPECT_EQ(tensor, expected);
}

// Rijndael's (AES) finite field GF(2^8)
// BuildMulTensor over F2[x]/(x^8+x^4+x^3+x+1).
TEST(BuildMulTensorTest, GF256) {
  // // x^8+x^4+x^3+x+1
  // p holds the coefficients of x^0..x^7; the x^8 term is implicit.
  std::bitset<8> p;
  p.set(0); // x^0
  p.set(1); // x^1
  p.set(3); // x^3
  p.set(4); // x^4
  const Tensor<2, 1, 8, 8, 8> tensor = BuildMulTensor<8>(p);
  // Example
  // https://en.wikipedia.org/wiki/Finite_field_arithmetic#Rijndael's_(AES)_finite_field
  // (x6 + x4 + x + 1) * (x7 + x6 + x3 + x) = 1
  std::array<uint8_t, 8> a = {1, 1, 0, 0, 1, 0, 1, 0};
  std::array<uint8_t, 8> b = {0, 1, 0, 1, 0, 0, 1, 1};
  std::array<uint8_t, 8> c = {1, 0, 0, 0, 0, 0, 0, 0};
  // Apply the tensor: c_k = sum_{i,j} T[i][j][k] * a_i * b_j  (mod 2 in F2).
  std::array<uint8_t, 8> actual_c{};
  for (std::size_t k = 0; k < 8; ++k) {
    int sum = 0;
    for (std::size_t i = 0; i < 8; ++i) {
      for (std::size_t j = 0; j < 8; ++j) {
        sum += tensor[i][j][k].value * a[i] * b[j];
      }
    }
    actual_c[k] = sum % 2;
  }
  EXPECT_EQ(actual_c, c);
}

// With p = x^n (all low-degree coefficients zero), x^n reduces to 0, so terms
// with i+j >= n vanish: ExtensionField equals Truncated multiplication.
TEST(BuildMulTensorTest, EqualsTruncatedWhenPIsXn) {
  const std::bitset<3> p3; // x^3
  EXPECT_EQ(BuildMulTensor<3>(p3), (truncated::BuildMulTensor<2, 1, 3>()));
  const std::bitset<5> p5; // x^5
  EXPECT_EQ(BuildMulTensor<5>(p5), (truncated::BuildMulTensor<2, 1, 5>()));
}

// With p = x^n - 1 = x^n + 1 (in F2), x^n reduces to 1, so degrees wrap around
// modulo n: ExtensionField equals Cyclic multiplication.
TEST(BuildMulTensorTest, EqualsCyclicWhenPIsXnMinus1) {
  std::bitset<3> p3;
  p3.set(0); // x^3 - 1
  EXPECT_EQ(BuildMulTensor<3>(p3), (cyclic::BuildMulTensor<2, 1, 3>()));
  std::bitset<5> p5;
  p5.set(0); // x^5 - 1
  EXPECT_EQ(BuildMulTensor<5>(p5), (cyclic::BuildMulTensor<2, 1, 5>()));
}

// --- Tests for the general (P, M) overload BuildMulTensor<P, M, N> ---

namespace {

// Apply the tensor as a bilinear product: c_k = Σ_{i,j} T[i][j][k] · a_i · b_j,
// all arithmetic in 𝔽_q.
template <int P, int M, std::size_t N>
std::array<GF<P, M>, N>
MulViaTensor(const Tensor<P, M, N, N, N> &t,
             const std::array<GF<P, M>, N> &a,
             const std::array<GF<P, M>, N> &b) {
  using F = GF<P, M>;
  std::array<F, N> c{};
  for (std::size_t i = 0; i < N; ++i) {
    if (a[i].value == 0)
      continue;
    for (std::size_t j = 0; j < N; ++j) {
      if (b[j].value == 0)
        continue;
      const F ab = F::Mul(a[i], b[j]);
      for (std::size_t k = 0; k < N; ++k) {
        c[k] = F::Add(c[k], F::Mul(t[i][j][k], ab));
      }
    }
  }
  return c;
}

// The idx-th element of 𝔽_q^N in the q-ary encoding (low coordinate first).
template <int P, int M, std::size_t N>
std::array<GF<P, M>, N> ElemAt(uint64_t idx) {
  const GFVec<P, M, static_cast<int>(N)> v =
      DecodeGFVec<P, M, static_cast<int>(N)>(idx);
  std::array<GF<P, M>, N> a{};
  for (std::size_t k = 0; k < N; ++k)
    a[k] = v[static_cast<int>(k)];
  return a;
}

// Verify the tensor built from the lex-smallest irreducible p makes 𝔽_q^N a
// field: the induced product is commutative, has 1 = x^0 as a two-sided
// identity, is associative, and every nonzero element is invertible. A finite
// commutative unital associative ring of q^N elements with no zero divisors is
// the field 𝔽_{q^N}, so this pins the tensor down independently of which
// irreducible was chosen (no hardcoded reduction table).
template <int P, int M, std::size_t N> void CheckTensorIsField() {
  using F = GF<P, M>;
  const auto p_low =
      extfield::IrreduciblePolyCoeffs<P, M, static_cast<int>(N)>();
  const Tensor<P, M, N, N, N> t = BuildMulTensor<P, M, N>(p_low);
  const uint64_t qn = IntPow(F::kQ, static_cast<int>(N));
  std::array<F, N> one{};
  one[0] = F::One();

  for (uint64_t ia = 0; ia < qn; ++ia) {
    const auto a = ElemAt<P, M, N>(ia);
    EXPECT_EQ((MulViaTensor<P, M, N>(t, one, a)), a) << "identity ia=" << ia;
    EXPECT_EQ((MulViaTensor<P, M, N>(t, a, one)), a) << "identity ia=" << ia;
    bool invertible = (ia == 0); // 0 needs no inverse
    for (uint64_t ib = 0; ib < qn; ++ib) {
      const auto b = ElemAt<P, M, N>(ib);
      const auto ab = MulViaTensor<P, M, N>(t, a, b);
      EXPECT_EQ(ab, (MulViaTensor<P, M, N>(t, b, a)))
          << "commutativity ia=" << ia << " ib=" << ib;
      if (ab == one)
        invertible = true;
    }
    EXPECT_TRUE(invertible) << "no inverse for element " << ia;
  }

  // Associativity over all triples (accumulated; q^N ≤ 81 in these tests).
  bool associative = true;
  for (uint64_t ia = 0; associative && ia < qn; ++ia) {
    const auto a = ElemAt<P, M, N>(ia);
    for (uint64_t ib = 0; associative && ib < qn; ++ib) {
      const auto b = ElemAt<P, M, N>(ib);
      const auto ab = MulViaTensor<P, M, N>(t, a, b);
      for (uint64_t ic = 0; ic < qn; ++ic) {
        const auto c = ElemAt<P, M, N>(ic);
        const auto lhs = MulViaTensor<P, M, N>(t, ab, c);
        const auto rhs =
            MulViaTensor<P, M, N>(t, a, MulViaTensor<P, M, N>(t, b, c));
        if (lhs != rhs) {
          associative = false;
          ADD_FAILURE() << "associativity ia=" << ia << " ib=" << ib
                        << " ic=" << ic;
          break;
        }
      }
    }
  }
  EXPECT_TRUE(associative);
}

// The GF<2, 1> overload must reproduce the trusted F_2 bitset path bit for bit.
template <std::size_t N> void CheckMatchesBitsetPath() {
  const auto bits = extfield::IrreduciblePolyBits<N>();
  std::array<GF<2, 1>, N> p_gf{};
  for (std::size_t i = 0; i < N; ++i) {
    p_gf[i] = GF<2, 1>{static_cast<uint8_t>(bits[i] ? 1 : 0)};
  }
  EXPECT_EQ((BuildMulTensor<N>(bits)),
            (BuildMulTensor<2, 1, N>(p_gf)))
      << "N=" << N;
}

} // namespace

TEST(BuildMulTensorFqTest, MatchesF2BitsetPathForM1) {
  CheckMatchesBitsetPath<2>();
  CheckMatchesBitsetPath<3>();
  CheckMatchesBitsetPath<4>();
  CheckMatchesBitsetPath<8>();
}

// M = 1 (prime base field) through the GF<P, M> overload, including odd P.
TEST(BuildMulTensorFqTest, DefinesFieldForM1) {
  CheckTensorIsField<2, 1, 2>(); // 𝔽_4
  CheckTensorIsField<2, 1, 3>(); // 𝔽_8
  CheckTensorIsField<3, 1, 2>(); // 𝔽_9
  CheckTensorIsField<5, 1, 2>(); // 𝔽_25
}

// M ≥ 2: the base field is itself an extension 𝔽_q.
TEST(BuildMulTensorFqTest, DefinesFieldGF16OverF4) {
  CheckTensorIsField<2, 2, 2>(); // 𝔽_{16} = 𝔽_{4^2}
}
TEST(BuildMulTensorFqTest, DefinesFieldGF64OverF4) {
  CheckTensorIsField<2, 2, 3>(); // 𝔽_{64} = 𝔽_{4^3}
}
TEST(BuildMulTensorFqTest, DefinesFieldGF81OverF9) {
  CheckTensorIsField<3, 2, 2>(); // 𝔽_{81} = 𝔽_{9^2}
}

} // namespace extfield
