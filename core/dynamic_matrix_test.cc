#include "dynamic_matrix.h"

#include <gtest/gtest.h>

TEST(DynamicMatrixTest, Constructor_ZeroMatrix) {
  DynamicMatrix<2, 1> m(2, 3);
  EXPECT_EQ(m.rows(), 2);
  EXPECT_EQ(m.cols(), 3);
  EXPECT_TRUE(m.IsZero());
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_EQ(m(i, j), 0);
    }
  }
  EXPECT_EQ(m.Rank(), 0);
}

TEST(DynamicMatrixTest, Resize) {
  DynamicMatrix<2, 1> m(2, 3);
  m(0, 0) = 1;
  m(1, 2) = 1;
  EXPECT_EQ(m.ToString(), "[1,0,0;0,0,1]");
  m.ResizeRows(3);
  EXPECT_EQ(m.rows(), 3);
  EXPECT_EQ(m.cols(), 3);
  EXPECT_EQ(m.ToString(), "[1,0,0;0,0,1;0,0,0]");
}

TEST(DynamicMatrixTest, SetAndGet) {
  DynamicMatrix<2, 1> m(2, 2);
  m(0, 0) = 1;
  m(1, 1) = 1;
  EXPECT_EQ(m(0, 0), 1);
  EXPECT_EQ(m(0, 1), 0);
  EXPECT_EQ(m(1, 0), 0);
  EXPECT_EQ(m(1, 1), 1);
  EXPECT_FALSE(m.IsZero());
  EXPECT_EQ(m.Rank(), 2);
}

TEST(DynamicMatrixTest, ToString) {
  DynamicMatrix<2, 1> m(3, 4);
  m(0, 0) = 1;
  m(0, 3) = 1; // row 0: 1001
  m(1, 0) = 1;
  m(1, 1) = 1;
  m(1, 3) = 1; // row 1: 1101
  m(2, 3) = 1; // row 2: 0001
  EXPECT_EQ(m.ToString(), "[1,0,0,1;1,1,0,1;0,0,0,1]");
}

TEST(DynamicMatrixTest, ToString_Zero) {
  DynamicMatrix<2, 1> m(2, 2);
  EXPECT_EQ(m.ToString(), "[0,0;0,0]");
}

TEST(DynamicMatrixTest, ToString_SingleFieldent) {
  DynamicMatrix<2, 1> m(1, 1);
  m(0, 0) = 1;
  EXPECT_EQ(m.ToString(), "[1]");
}

TEST(DynamicMatrixTest, Plus_SameSize) {
  DynamicMatrix<2, 1> a(2, 2);
  a(0, 0) = 1;
  a(1, 1) = 1;
  DynamicMatrix<2, 1> b(2, 2);
  b(0, 0) = 1;

  DynamicMatrix<2, 1> sum = a.Plus(b);
  EXPECT_EQ(sum(0, 0), 0); // 1 + 1 = 0 in F_2
  EXPECT_EQ(sum(1, 1), 1); // 0 + 1 = 1
  EXPECT_EQ(sum.rows(), 2);
  EXPECT_EQ(sum.cols(), 2);
}

TEST(DynamicMatrixTest, Plus_Zero) {
  DynamicMatrix<2, 1> a(2, 2);
  a(0, 1) = 1;
  DynamicMatrix<2, 1> b(2, 2);

  DynamicMatrix<2, 1> sum = a.Plus(b);
  EXPECT_EQ(sum(0, 1), 1);
  EXPECT_TRUE(sum.Plus(a).IsZero()); // sum + a = a + a = 0 in F_2
}

TEST(DynamicMatrixTest, IsZero) {
  DynamicMatrix<2, 1> m(2, 2);
  EXPECT_TRUE(m.IsZero());
  m(0, 0) = 1;
  EXPECT_FALSE(m.IsZero());
  m(0, 0) = 0;
  EXPECT_TRUE(m.IsZero());
}

TEST(DynamicMatrixTest, Rank_Zero) {
  DynamicMatrix<2, 1> m(3, 4);
  EXPECT_EQ(m.Rank(), 0);
}

TEST(DynamicMatrixTest, Rank_Identity2x2) {
  DynamicMatrix<2, 1> m(2, 2);
  m(0, 0) = 1;
  m(1, 1) = 1;
  EXPECT_EQ(m.Rank(), 2);
}

TEST(DynamicMatrixTest, Rank_Identity3x3) {
  DynamicMatrix<2, 1> m(3, 3);
  m(0, 0) = 1;
  m(1, 1) = 1;
  m(2, 2) = 1;
  EXPECT_EQ(m.Rank(), 3);
}

TEST(DynamicMatrixTest, Rank_One) {
  DynamicMatrix<2, 1> m(2, 2);
  m(0, 0) = 1;
  m(0, 1) = 1;
  // Row 0: [1, 1], Row 1: [0, 0]
  EXPECT_EQ(m.Rank(), 1);
}

TEST(DynamicMatrixTest, Rank_LinearlyDependentRows) {
  DynamicMatrix<2, 1> m(2, 2);
  m(0, 0) = 1;
  m(0, 1) = 1;
  m(1, 0) = 1;
  m(1, 1) = 1;
  // Both rows [1, 1]
  EXPECT_EQ(m.Rank(), 1);
}

TEST(DynamicMatrixTest, Rank_NonSquare_RowsMoreThanCols) {
  DynamicMatrix<2, 1> m(4, 2);
  m(0, 0) = 1;
  m(1, 1) = 1;
  // Two independent columns, rank = 2
  EXPECT_EQ(m.Rank(), 2);
}

TEST(DynamicMatrixTest, Rank_NonSquare_ColsMoreThanRows) {
  DynamicMatrix<2, 1> m(2, 4);
  m(0, 0) = 1;
  m(0, 1) = 1;
  m(1, 2) = 1;
  m(1, 3) = 1;
  // Two independent rows, rank = 2
  EXPECT_EQ(m.Rank(), 2);
}

TEST(DynamicMatrixTest, Rank_EmptyRows) {
  DynamicMatrix<2, 1> m(0, 3);
  EXPECT_EQ(m.Rank(), 0);
}

TEST(DynamicMatrixTest, Rank_EmptyCols) {
  DynamicMatrix<2, 1> m(3, 0);
  EXPECT_EQ(m.Rank(), 0);
}

TEST(DynamicMatrixTest, Plus_DimensionsMatch) {
  DynamicMatrix<2, 1> a(1, 1);
  DynamicMatrix<2, 1> b(1, 1);
  a(0, 0) = 1;
  DynamicMatrix<2, 1> sum = a.Plus(b);
  EXPECT_EQ(sum.ToString(), "[1]");
}

// --- F_q (M >= 2) coverage ---------------------------------------------------
//
// 𝔽_4 = 𝔽_2[y]/(y^2 + y + 1). Fieldents are F_q indices in [0, 4); the index
// decomposes into base-2 digits (a_0, a_1) low-coef first, so:
//   0 = 0 (digits {0,0}), 1 = 1 (digits {1,0}), y = 2 (digits {0,1}),
//   y+1 = 3 (digits {1,1}).
// y is the primitive: y · y = y+1, y · (y+1) = 1.

TEST(DynamicMatrixFqTest, Q4_RankPivotsViaInverse) {
  // Diagonal of two distinct non-1 nonzero entries: pivot normalization
  // via Inverse must put both rows in canonical RREF and report rank 2.
  using M = DynamicMatrix<2, 2>;
  using Field = M::Field;
  const Field y = 2;
  const Field y_plus_1 = 3;
  M mat(2, 2);
  mat(0, 0) = y;
  mat(1, 1) = y_plus_1;
  EXPECT_EQ(mat.Rank(), 2);
}

TEST(DynamicMatrixFqTest, Q4_RankDetectsDependentRow) {
  // Row 1 = y · Row 0 over 𝔽_4. y · y = y+1, y · 1 = y.
  using M = DynamicMatrix<2, 2>;
  using Field = M::Field;
  const Field one = 1;
  const Field y = 2;
  const Field y_plus_1 = 3;
  M mat(2, 2);
  mat(0, 0) = y;
  mat(0, 1) = one;
  mat(1, 0) = y_plus_1; // = y · y
  mat(1, 1) = y;        // = y · 1
  EXPECT_EQ(mat.Rank(), 1);
}

TEST(DynamicMatrixFqTest, Q4_PlusIsFieldAddition) {
  using M = DynamicMatrix<2, 2>;
  using Field = M::Field;
  M a(1, 2);
  M b(1, 2);
  a(0, 0) = Field{2}; // y
  a(0, 1) = Field{3}; // y+1
  b(0, 0) = Field{3}; // y+1
  b(0, 1) = Field{3}; // y+1
  M sum = a.Plus(b);
  EXPECT_EQ(sum(0, 0), Field{1}); // y + (y+1) = 1
  EXPECT_EQ(sum(0, 1), Field{0}); // (y+1) + (y+1) = 0
}

TEST(DynamicMatrixFqTest, Q4_ToStringDigitOrder) {
  using M = DynamicMatrix<2, 2>;
  using Field = M::Field;
  M mat(2, 3);
  mat(0, 0) = Field{1}; // 1   -> "10"
  mat(0, 1) = Field{2}; // y   -> "01"
  mat(0, 2) = Field{3}; // y+1 -> "11"
  mat(1, 0) = Field{1}; // 1   -> "10"
  mat(1, 1) = Field{2}; // y   -> "01"
  mat(1, 2) = Field{0}; // 0   -> "00"
  EXPECT_EQ(mat.ToString(), "[10,01,11;10,01,00]");
}
