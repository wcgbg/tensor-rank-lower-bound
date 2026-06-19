#include "problems/cyclic/tensor.h"

#include <string>

#include <gtest/gtest.h>

#include "core/tensor.h"
#include "core/tensor_utils.h"

// BuildMulTensor: T[i][j][(i+j) % n] = 1.
TEST(BuildMulTensorTest, N2) {
  // (i+j) mod 2: a1*b1 wraps back to c0.
  EXPECT_EQ(
      (TensorToSparseString<2, 1, 2, 2, 2>(cyclic::BuildMulTensor<2, 1, 2>())),
      "a0*b0*c0 + a0*b1*c1 + a1*b0*c1 + a1*b1*c0");
}

TEST(BuildMulTensorTest, N3) {
  using T = Tensor<2, 1, 3, 3, 3>;
  const T tensor = cyclic::BuildMulTensor<2, 1, 3>();
  std::string expected_string;
  expected_string += "+a0*b0*c0+a1*b2*c0+a2*b1*c0";
  expected_string += "+a0*b1*c1+a1*b0*c1+a2*b2*c1";
  expected_string += "+a0*b2*c2+a1*b1*c2+a2*b0*c2";
  EXPECT_EQ((SparseStringToTensor<2, 1, 3, 3, 3>(expected_string)), tensor);
}
