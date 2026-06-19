#include "problems/full/tensor.h"

#include <string>

#include <gtest/gtest.h>

#include "core/tensor.h"
#include "core/tensor_utils.h"

// BuildMulTensor: T[i][j][i+j] = 1.
TEST(BuildMulTensorTest, N2) {
  using T = Tensor<2, 1, 2, 2, 3>;
  const T tensor = full::BuildMulTensor<2, 1, 2>();
  EXPECT_EQ((TensorToSparseString<2, 1, 2, 2, 3>(tensor)),
            "a0*b0*c0 + a0*b1*c1 + a1*b0*c1 + a1*b1*c2");
}

TEST(BuildMulTensorTest, N3) {
  using T = Tensor<2, 1, 3, 3, 5>;
  const T actual_tensor = full::BuildMulTensor<2, 1, 3>();
  std::string expected_string;
  expected_string += "+a0*b0*c0";
  expected_string += "+a0*b1*c1+a1*b0*c1";
  expected_string += "+a0*b2*c2+a1*b1*c2+a2*b0*c2";
  expected_string += "+a1*b2*c3+a2*b1*c3";
  expected_string += "+a2*b2*c4";
  EXPECT_EQ((SparseStringToTensor<2, 1, 3, 3, 5>(expected_string)),
            actual_tensor);
}
