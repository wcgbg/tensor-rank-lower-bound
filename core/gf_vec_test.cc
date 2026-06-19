#include "core/gf_vec.h"

#include <gtest/gtest.h>

namespace {

TEST(GFVecF3Test, ConstructionAndAccess) {
  GFVec<3, 1, 4> v{};
  EXPECT_TRUE(v.IsZero());
  v[0] = 1;
  v[3] = 2;
  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[3], 2);
  EXPECT_FALSE(v.IsZero());
}

TEST(GFVecF3Test, AdditionWraps) {
  GFVec<3, 1, 3> a{{1, 2, 0}};
  GFVec<3, 1, 3> b{{2, 2, 1}};
  const GFVec<3, 1, 3> sum = a + b;
  EXPECT_EQ(sum[0], 0); // 1 + 2 = 0 mod 3
  EXPECT_EQ(sum[1], 1); // 2 + 2 = 1 mod 3
  EXPECT_EQ(sum[2], 1);
}

TEST(GFVecF3Test, SubtractionWraps) {
  GFVec<3, 1, 3> a{{0, 0, 0}};
  GFVec<3, 1, 3> b{{1, 2, 0}};
  const GFVec<3, 1, 3> diff = a - b;
  EXPECT_EQ(diff[0], 2); // 0 - 1 = 2 mod 3
  EXPECT_EQ(diff[1], 1); // 0 - 2 = 1 mod 3
  EXPECT_EQ(diff[2], 0);
}

TEST(GFVecF3Test, NegationIsAdditiveInverse) {
  GFVec<3, 1, 3> a{{1, 2, 0}};
  const GFVec<3, 1, 3> sum = a + (-a);
  EXPECT_TRUE(sum.IsZero());
}

TEST(GFVecF3Test, ScalarMul) {
  GFVec<3, 1, 3> a{{1, 2, 0}};
  const GFVec<3, 1, 3> doubled = uint8_t{2} * a;
  EXPECT_EQ(doubled[0], 2);
  EXPECT_EQ(doubled[1], 1); // 2*2 = 4 = 1 mod 3
  EXPECT_EQ(doubled[2], 0);
}

TEST(GFVecF3Test, LeadingNonzeroIdx) {
  GFVec<3, 1, 5> v{};
  EXPECT_EQ(v.LeadingNonzeroIdx(), -1);
  v[2] = 1;
  EXPECT_EQ(v.LeadingNonzeroIdx(), 2);
  v[4] = 2;
  EXPECT_EQ(v.LeadingNonzeroIdx(), 4);
}

TEST(GFVecF3Test, Equality) {
  GFVec<3, 1, 3> a{{1, 2, 0}};
  GFVec<3, 1, 3> b{{1, 2, 0}};
  GFVec<3, 1, 3> c{{1, 2, 1}};
  EXPECT_EQ(a, b);
  EXPECT_NE(a, c);
}

TEST(GFVecF5Test, ConstexprOperations) {
  // Confirm the constexpr-ness of the value-semantics operators.
  static constexpr GFVec<5, 1, 3> a{{3, 4, 2}};
  static constexpr GFVec<5, 1, 3> b{{4, 1, 3}};
  static constexpr GFVec<5, 1, 3> sum = a + b;
  static_assert(sum[0] == 2); // 3 + 4 = 7 = 2 mod 5
  static_assert(sum[1] == 0); // 4 + 1 = 5 = 0 mod 5
  static_assert(sum[2] == 0); // 2 + 3 = 5 = 0 mod 5
  SUCCEED();
}

TEST(EncodeGFVecTest, MatchesBaseExpansion) {
  // F_3, N=3, M=1: v = (a_0, a_1, a_2) ↦ a_0 + 3 a_1 + 9 a_2.
  EXPECT_EQ((EncodeGFVec<3, 1, 3>(GFVec<3, 1, 3>{{0, 0, 0}})), 0u);
  EXPECT_EQ((EncodeGFVec<3, 1, 3>(GFVec<3, 1, 3>{{1, 0, 0}})), 1u);
  EXPECT_EQ((EncodeGFVec<3, 1, 3>(GFVec<3, 1, 3>{{2, 0, 0}})), 2u);
  EXPECT_EQ((EncodeGFVec<3, 1, 3>(GFVec<3, 1, 3>{{0, 1, 0}})), 3u);
  EXPECT_EQ((EncodeGFVec<3, 1, 3>(GFVec<3, 1, 3>{{1, 1, 0}})), 4u);
  EXPECT_EQ((EncodeGFVec<3, 1, 3>(GFVec<3, 1, 3>{{2, 2, 2}})), 26u);
}

TEST(EncodeGFVecTest, RoundTripsWithDecode) {
  // Every integer in [0, P^N) round-trips through Decode∘Encode.
  constexpr uint64_t kRange = 3 * 3 * 3 * 3; // 3^4
  for (uint64_t i = 0; i < kRange; ++i) {
    const auto v = DecodeGFVec<3, 1, 4>(i);
    EXPECT_EQ((EncodeGFVec<3, 1, 4>(v)), i) << "i=" << i;
  }
  // And every GFVec round-trips through Encode∘Decode.
  for (uint8_t a = 0; a < 5; ++a) {
    for (uint8_t b = 0; b < 5; ++b) {
      const GFVec<5, 1, 2> v{{a, b}};
      EXPECT_EQ((DecodeGFVec<5, 1, 2>(EncodeGFVec<5, 1, 2>(v))), v)
          << "a=" << int{a} << " b=" << int{b};
    }
  }
}

TEST(EncodeGFVecTest, OrderingMatchesIntegerCompare) {
  // The operator<=> definition promises that integer order ↔ GFVec lex order
  // when the integer is the EncodeGFVec value. Spot-check across F_3, N=3.
  for (uint64_t i = 0; i < 27; ++i) {
    for (uint64_t j = 0; j < 27; ++j) {
      const auto vi = DecodeGFVec<3, 1, 3>(i);
      const auto vj = DecodeGFVec<3, 1, 3>(j);
      EXPECT_EQ(i < j, vi < vj) << "i=" << i << " j=" << j;
      const bool lt = (EncodeGFVec<3, 1, 3>(vi) < EncodeGFVec<3, 1, 3>(vj));
      EXPECT_EQ(lt, vi < vj);
    }
  }
}

TEST(EncodeGFVecTest, ConstexprEncode) {
  static constexpr GFVec<5, 1, 3> v{{4, 0, 3}};
  static_assert(EncodeGFVec<5, 1, 3>(v) == 4u + 0u * 5 + 3u * 25);
  SUCCEED();
}

// --- F_q (M >= 2) layout tests ---

TEST(GFVecF4Test, LayoutAndAccess) {
  // F_4 = F_2[y]/g(y) with g(y) = y^2 + y + 1; M = 2, q = 4.
  // GFVec<2, 2, 3> now stores 3 F_q indices = 3 bytes (one per element).
  static_assert(sizeof(GFVec<2, 2, 3>) == 3);
  GFVec<2, 2, 3> v{};
  EXPECT_TRUE(v.IsZero());
  // Set element 0 to "y + 1" = the primitive of F_4*. Digit-vector (1, 1)
  // encodes to the index 1 + 1·2 = 3.
  v.Set(0, 3);
  EXPECT_FALSE(v.IsZero());
  EXPECT_EQ(v[0], 3);
  EXPECT_EQ(v[1], 0);
  EXPECT_EQ(v.LeadingNonzeroIdx(), 0);
}

TEST(EncodeGFVecTest, RoundTripsM2) {
  // q=4, N=3: 4^3 = 64 values.
  for (uint64_t i = 0; i < 64; ++i) {
    const auto v = DecodeGFVec<2, 2, 3>(i);
    EXPECT_EQ((EncodeGFVec<2, 2, 3>(v)), i) << "i=" << i;
  }
}

} // namespace
