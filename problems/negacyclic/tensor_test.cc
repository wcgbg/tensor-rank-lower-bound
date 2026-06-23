#include "problems/negacyclic/tensor.h"

#include <gtest/gtest.h>

#include "core/gf.h"
#include "core/tensor.h"
#include "core/tensor_utils.h"
#include "problems/cyclic/tensor.h"

namespace negacyclic {
namespace {

// Over characteristic 2, x^N + 1 = x^N - 1, so the negacyclic tensor is
// identical to the cyclic tensor.
TEST(BuildMulTensorTest, EqualsCyclicOverF2) {
  EXPECT_EQ((BuildMulTensor<2, 1, 2>()), (cyclic::BuildMulTensor<2, 1, 2>()));
  EXPECT_EQ((BuildMulTensor<2, 1, 3>()), (cyclic::BuildMulTensor<2, 1, 3>()));
  EXPECT_EQ((BuildMulTensor<2, 1, 5>()), (cyclic::BuildMulTensor<2, 1, 5>()));
  EXPECT_EQ((BuildMulTensor<2, 4, 4>()), (cyclic::BuildMulTensor<2, 4, 4>()));
}

// Over F_3, x^2 + 1: x^2 ≡ -1, so a1*b1 lands at c0 with coefficient -1 = 2.
TEST(BuildMulTensorTest, F3N2Signs) {
  using F = GF<3, 1>;
  const auto t = BuildMulTensor<3, 1, 2>();
  const F one = F::One();
  const F neg_one = F{static_cast<uint8_t>(2)}; // -1 mod 3

  // a0*b0 -> c0 (+1); a0*b1 -> c1 (+1); a1*b0 -> c1 (+1); a1*b1 -> c0 (-1).
  EXPECT_EQ(t[0][0][0], one);
  EXPECT_EQ(t[0][0][1], F::Zero());
  EXPECT_EQ(t[0][1][1], one);
  EXPECT_EQ(t[1][0][1], one);
  EXPECT_EQ(t[1][1][0], neg_one);
  EXPECT_EQ(t[1][1][1], F::Zero());

  EXPECT_EQ(TensorToSparseString(t),
            "a0*b0*c0 + a0*b1*c1 + a1*b0*c1 + 2*a1*b1*c0");
}

// Over F_5, x^3 + 1: x^3 ≡ -1 = 4. a2*b1 (i+j=3) -> c0 with coefficient -1.
TEST(BuildMulTensorTest, F5N3Wraparound) {
  using F = GF<5, 1>;
  const auto t = BuildMulTensor<5, 1, 3>();
  const F neg_one = F{static_cast<uint8_t>(4)}; // -1 mod 5

  EXPECT_EQ((t[0][0][0]), F::One()); // i+j=0 < 3
  EXPECT_EQ((t[1][1][2]), F::One()); // i+j=2 < 3
  EXPECT_EQ((t[2][1][0]), neg_one);  // i+j=3 >= 3 -> c0, -1
  EXPECT_EQ((t[2][2][1]), neg_one);  // i+j=4 >= 3 -> c1, -1

  std::string expected_string;
  expected_string += "+ a0*b0*c0 + 4*a1*b2*c0 + 4*a2*b1*c0";
  expected_string += "+ a0*b1*c1 + a1*b0*c1 + 4*a2*b2*c1";
  expected_string += "+ a0*b2*c2 + a1*b1*c2 + a2*b0*c2";
  EXPECT_EQ((SparseStringToTensor<5, 1, 3, 3, 3>(expected_string)), t);
}

} // namespace
} // namespace negacyclic
