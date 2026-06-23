#include "problems/truncated/tensor.h"

#include <string>

#include <gtest/gtest.h>

#include "core/tensor.h"
#include "core/tensor_utils.h"

// BuildMulTensor: T[i][j][i+j] = 1 iff i+j < n.
TEST(BuildMulTensorTest, N2) {
  // a1*b1 would land at c2, which is truncated.
  using T = Tensor<2, 1, 2, 2, 2>;
  const T tensor = truncated::BuildMulTensor<2, 1, 2>();
  EXPECT_EQ((TensorToSparseString<2, 1, 2, 2, 2>(tensor)),
            "a0*b0*c0 + a0*b1*c1 + a1*b0*c1");
}

TEST(BuildMulTensorTest, N3) {
  using T = Tensor<2, 1, 3, 3, 3>;
  const T tensor = truncated::BuildMulTensor<2, 1, 3>();
  std::string expected_string;
  expected_string += "+a0*b0*c0";
  expected_string += "+a0*b1*c1+a1*b0*c1";
  expected_string += "+a0*b2*c2+a1*b1*c2+a2*b0*c2";
  EXPECT_EQ((SparseStringToTensor<2, 1, 3, 3, 3>(expected_string)), tensor);
}
