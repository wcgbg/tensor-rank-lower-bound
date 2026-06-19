#include "tensor_isomorphism.h"

#include <cstddef>
#include <cstdint>
#include <random>

#include <gtest/gtest.h>

#include "core/gf.h"
#include "core/gf_vec.h"
#include "tensor_utils.h"

namespace {

using tensor_isomorphism_internal::ApplyToModeA;
using tensor_isomorphism_internal::ApplyToModeB;
using tensor_isomorphism_internal::ApplyToModeC;
using tensor_isomorphism_internal::EnumerateInvertibleMatrices;
using tensor_isomorphism_internal::SquareMatrix;

// Brute-force GL(NA) × GL(NB) × GL(NC) oracle — the original implementation of
// AreTensorsIsomorphic, kept here as an independent ground truth for the
// randomized cross-check tests below. Viable only at the tiny shapes the tests
// exercise.
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
bool AreTensorsIsomorphicSlow(const Tensor<P, M, NA, NB, NC> &t0,
                              const Tensor<P, M, NA, NB, NC> &t1) {
  const auto group_a = EnumerateInvertibleMatrices<P, M, NA>();
  const auto group_b = EnumerateInvertibleMatrices<P, M, NB>();
  const auto group_c = EnumerateInvertibleMatrices<P, M, NC>();
  for (const auto &a : group_a) {
    const Tensor<P, M, NA, NB, NC> ta = ApplyToModeA(t0, a);
    for (const auto &b : group_b) {
      const Tensor<P, M, NA, NB, NC> tab = ApplyToModeB(ta, b);
      for (const auto &c : group_c) {
        if (ApplyToModeC(tab, c) == t1) {
          return true;
        }
      }
    }
  }
  return false;
}

template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
Tensor<P, M, NA, NB, NC> RandomTensor(std::mt19937_64 &rng) {
  Tensor<P, M, NA, NB, NC> t{};
  std::uniform_int_distribution<int> dist(0, GF<P, M>::kQ - 1);
  for (std::size_t i = 0; i < NA; ++i) {
    for (std::size_t j = 0; j < NB; ++j) {
      for (std::size_t k = 0; k < NC; ++k) {
        t[i][j][k] = GF<P, M>{static_cast<uint8_t>(dist(rng))};
      }
    }
  }
  return t;
}

template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
void CrossCheckRandom(int num_trials, std::mt19937_64 &rng) {
  const auto group_a = EnumerateInvertibleMatrices<P, M, NA>();
  const auto group_b = EnumerateInvertibleMatrices<P, M, NB>();
  const auto group_c = EnumerateInvertibleMatrices<P, M, NC>();
  std::uniform_int_distribution<std::size_t> pa(0, group_a.size() - 1);
  std::uniform_int_distribution<std::size_t> pb(0, group_b.size() - 1);
  std::uniform_int_distribution<std::size_t> pc(0, group_c.size() - 1);
  for (int trial = 0; trial < num_trials; ++trial) {
    const auto t0 = RandomTensor<P, M, NA, NB, NC>(rng);
    const auto t1 = RandomTensor<P, M, NA, NB, NC>(rng);
    EXPECT_EQ(AreTensorsIsomorphic(t0, t1), AreTensorsIsomorphicSlow(t0, t1))
        << "P=" << P << " M=" << M << " shape=(" << NA << "," << NB << "," << NC
        << ") trial=" << trial << " (independent random pair)";

    const auto &a = group_a[pa(rng)];
    const auto &b = group_b[pb(rng)];
    const auto &c = group_c[pc(rng)];
    const auto t1_iso = ApplyToModeC(ApplyToModeB(ApplyToModeA(t0, a), b), c);
    EXPECT_TRUE(AreTensorsIsomorphic(t0, t1_iso))
        << "P=" << P << " M=" << M << " shape=(" << NA << "," << NB << "," << NC
        << ") trial=" << trial << " (constructed isomorphic pair, fast)";
    EXPECT_TRUE(AreTensorsIsomorphicSlow(t0, t1_iso))
        << "P=" << P << " M=" << M << " shape=(" << NA << "," << NB << "," << NC
        << ") trial=" << trial << " (constructed isomorphic pair, slow)";
  }
}

// |GL(n, F_2)| = prod_{i=0}^{n-1} (2^n - 2^i).
TEST(EnumerateInvertibleMatricesTest, GroupOrdersF2) {
  EXPECT_EQ((EnumerateInvertibleMatrices<2, 1, 1>().size()), 1u);
  EXPECT_EQ((EnumerateInvertibleMatrices<2, 1, 2>().size()), 6u);
  EXPECT_EQ((EnumerateInvertibleMatrices<2, 1, 3>().size()), 168u);
  EXPECT_EQ((EnumerateInvertibleMatrices<2, 1, 4>().size()), 20160u);
}

// |GL(n, F_3)| = prod_{i=0}^{n-1} (3^n - 3^i).
TEST(EnumerateInvertibleMatricesTest, GroupOrdersF3) {
  EXPECT_EQ((EnumerateInvertibleMatrices<3, 1, 1>().size()), 2u);
  EXPECT_EQ((EnumerateInvertibleMatrices<3, 1, 2>().size()), 48u);
  EXPECT_EQ((EnumerateInvertibleMatrices<3, 1, 3>().size()), 11232u);
}

// |GL(n, F_4)| = prod_{i=0}^{n-1} (4^n - 4^i). Forces the GFVec primary
// template and the F_4 Cayley tables (no BitVec specialization).
TEST(EnumerateInvertibleMatricesTest, GroupOrdersF4) {
  EXPECT_EQ((EnumerateInvertibleMatrices<2, 2, 1>().size()), 3u);
  EXPECT_EQ((EnumerateInvertibleMatrices<2, 2, 2>().size()), 180u);
}

TEST(AreTensorsIsomorphicTest, IdentityIsomorphism) {
  // Every tensor is isomorphic to itself (take A = B = C = identity).
  Tensor<2, 1, 3, 3, 3> t = SparseStringToTensor<2, 1, 3, 3, 3>(
      "a0*b0*c0 + a1*b2*c1 + a2*b1*c2 + a0*b1*c1");
  EXPECT_TRUE(AreTensorsIsomorphic(t, t));
}

TEST(AreTensorsIsomorphicTest, ZeroAndNonzeroAreNotIsomorphic) {
  // A linear map sends 0 to 0, so the zero tensor is alone in its class.
  Tensor<2, 1, 2, 2, 2> zero{};
  Tensor<2, 1, 2, 2, 2> nonzero{};
  nonzero[0][0][0] = 1;
  EXPECT_FALSE(AreTensorsIsomorphic(zero, nonzero));
  EXPECT_TRUE(AreTensorsIsomorphic(zero, zero));
}

TEST(AreTensorsIsomorphicTest, RecognizesAModeChangeOfBasis) {
  // Build t2 from t1 by acting with a single invertible matrix on mode A; the
  // two must come out isomorphic.
  Tensor<2, 1, 2, 2, 2> t1{};
  t1[0][0][0] = 1;
  t1[1][1][1] = 1;
  // A = [[1,1],[0,1]] over F_2 (rows are bit patterns 0b11 and 0b10).
  SquareMatrix<2, 1, 2> a = {{{0b11}, {0b10}}};
  const Tensor<2, 1, 2, 2, 2> t2 = ApplyToModeA(t1, a);
  EXPECT_NE(t1, t2);
  EXPECT_TRUE(AreTensorsIsomorphic(t1, t2));
}

TEST(AreTensorsIsomorphicTest, RecognizesAllThreeModeChangeOfBasis) {
  Tensor<2, 1, 2, 2, 3> t1{};
  t1[0][1][2] = 1;
  t1[1][0][0] = 1;
  SquareMatrix<2, 1, 2> a = {{{0b11}, {0b01}}};
  SquareMatrix<2, 1, 2> b = {{{0b10}, {0b11}}};
  SquareMatrix<2, 1, 3> c = {{{0b011}, {0b010}, {0b101}}};
  const Tensor<2, 1, 2, 2, 3> t2 =
      ApplyToModeC(ApplyToModeB(ApplyToModeA(t1, a), b), c);
  EXPECT_TRUE(AreTensorsIsomorphic(t1, t2));
}

TEST(AreTensorsIsomorphicTest, DifferentRankNotIsomorphic) {
  // Tensor rank is an isomorphism invariant: a rank-1 tensor a0*b0*c0 cannot be
  // isomorphic to the rank-2 tensor a0*b0*c0 + a1*b1*c1.
  Tensor<2, 1, 2, 2, 2> rank1{};
  rank1[0][0][0] = 1;
  Tensor<2, 1, 2, 2, 2> rank2{};
  rank2[0][0][0] = 1;
  rank2[1][1][1] = 1;
  EXPECT_FALSE(AreTensorsIsomorphic(rank1, rank2));
}

TEST(AreTensorsIsomorphicTest, IsomorphismIsSymmetric) {
  Tensor<2, 1, 2, 2, 2> t1{};
  t1[0][0][0] = 1;
  t1[0][1][1] = 1;
  SquareMatrix<2, 1, 2> b = {{{0b11}, {0b10}}};
  const Tensor<2, 1, 2, 2, 2> t2 = ApplyToModeB(t1, b);
  EXPECT_TRUE(AreTensorsIsomorphic(t1, t2));
  EXPECT_TRUE(AreTensorsIsomorphic(t2, t1));
}

TEST(AreTensorsIsomorphicTest, F3ZeroAndNonzeroAreNotIsomorphic) {
  Tensor<3, 1, 2, 2, 2> zero{};
  Tensor<3, 1, 2, 2, 2> nonzero{};
  nonzero[0][0][0] = 2;
  EXPECT_FALSE(AreTensorsIsomorphic(zero, nonzero));
  EXPECT_TRUE(AreTensorsIsomorphic(zero, zero));
}

TEST(AreTensorsIsomorphicTest, F3RecognizesAModeChangeOfBasis) {
  // Two nonzero cells, then act by an invertible 2x2 over F_3 on mode A.
  Tensor<3, 1, 2, 2, 2> t1{};
  t1[0][0][0] = 1;
  t1[1][1][1] = 2;
  // A = [[2, 1], [0, 1]] over F_3: det = 2 != 0, so invertible.
  using V = GFVec<3, 1, 2>;
  SquareMatrix<3, 1, 2> a = {
      {V{{GF<3, 1>{2}, GF<3, 1>{1}}}, V{{GF<3, 1>{0}, GF<3, 1>{1}}}}};
  const Tensor<3, 1, 2, 2, 2> t2 = ApplyToModeA(t1, a);
  EXPECT_NE(t1, t2);
  EXPECT_TRUE(AreTensorsIsomorphic(t1, t2));
}

TEST(AreTensorsIsomorphicTest, F4RecognizesBModeChangeOfBasis) {
  // F_4 = GF<2, 2>: elements 0, 1, y (=2), y+1 (=3).
  // Single nonzero cell with value y; act by an invertible 2x2 on mode B.
  Tensor<2, 2, 2, 2, 2> t1{};
  t1[0][0][0] = 2; // y
  // B = [[1, y], [y, 1]]: det = 1 - y² = 1 - (y+1) = -y = y in char 2,
  // which is nonzero, so B is invertible.
  using V = GFVec<2, 2, 2>;
  using F = GF<2, 2>;
  SquareMatrix<2, 2, 2> b = {{V{{F{1}, F{2}}}, V{{F{2}, F{1}}}}};
  const Tensor<2, 2, 2, 2, 2> t2 = ApplyToModeB(t1, b);
  EXPECT_NE(t1, t2);
  EXPECT_TRUE(AreTensorsIsomorphic(t1, t2));
}

// Randomized cross-check: at each shape, draw random tensors and assert the
// fast AreTensorsIsomorphic agrees with the brute-force
// AreTensorsIsomorphicSlow oracle. Each trial also constructs a deliberately
// isomorphic pair (by applying a random (A, B, C) triple) and asserts both
// functions return true. Trial counts are tuned so each TEST runs in well under
// a second under
// --config=opt; the slow oracle's cost is
// |GL(NA,𝔽_q)|·|GL(NB,𝔽_q)|·|GL(NC,𝔽_q)| matrix-applies per call. The mt19937
// seed is fixed so failures reproduce.
TEST(AreTensorsIsomorphicTest, RandomCrossCheckF2_222) {
  std::mt19937_64 rng;
  CrossCheckRandom<2, 1, 2, 2, 2>(50, rng); // 6^3 = 216 triples / slow call
}
TEST(AreTensorsIsomorphicTest, RandomCrossCheckF2_223) {
  std::mt19937_64 rng;
  CrossCheckRandom<2, 1, 2, 2, 3>(20, rng); // 6·6·168 ≈ 6k triples / slow call
}
TEST(AreTensorsIsomorphicTest, RandomCrossCheckF2_333) {
  std::mt19937_64 rng;
  CrossCheckRandom<2, 1, 3, 3, 3>(3, rng); // 168^3 ≈ 4.7M triples / slow call
}
TEST(AreTensorsIsomorphicTest, RandomCrossCheckF3_222) {
  std::mt19937_64 rng;
  CrossCheckRandom<3, 1, 2, 2, 2>(10, rng); // 48^3 ≈ 110k triples / slow call
}
TEST(AreTensorsIsomorphicTest, RandomCrossCheckF4_222) {
  std::mt19937_64 rng;
  CrossCheckRandom<2, 2, 2, 2, 2>(3, rng); // 180^3 ≈ 5.8M triples / slow call
}

} // namespace
