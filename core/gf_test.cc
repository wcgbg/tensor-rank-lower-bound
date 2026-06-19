#include "core/gf.h"

#include <algorithm>

#include <gtest/gtest.h>

namespace {

TEST(GF2Test, AddSubXor) {
  using F = GF<2, 1>;
  EXPECT_EQ(F::Add(F{0}, F{0}), F{0});
  EXPECT_EQ(F::Add(F{1}, F{0}), F{1});
  EXPECT_EQ(F::Add(F{1}, F{1}), F{0});
  EXPECT_EQ(F::Sub(F{0}, F{1}), F{1});
  EXPECT_EQ(F::Sub(F{1}, F{1}), F{0});
}

TEST(GF2Test, MulNeg) {
  using F = GF<2, 1>;
  EXPECT_EQ(F::Mul(F{0}, F{0}), F{0});
  EXPECT_EQ(F::Mul(F{1}, F{1}), F{1});
  EXPECT_EQ(F::Mul(F{1}, F{0}), F{0});
  EXPECT_EQ(F::Neg(F{0}), F{0});
  EXPECT_EQ(F::Neg(F{1}), F{1});
}

TEST(GF2Test, Inverse) {
  using F = GF<2, 1>;
  EXPECT_EQ(F::Inverse(F{1}), F{1});
}

TEST(GF3Test, AddSub) {
  using F = GF<3, 1>;
  EXPECT_EQ(F::Add(F{1}, F{2}), F{0});
  EXPECT_EQ(F::Add(F{2}, F{2}), F{1});
  EXPECT_EQ(F::Sub(F{0}, F{1}), F{2});
  EXPECT_EQ(F::Sub(F{1}, F{2}), F{2});
}

TEST(GF3Test, MulNeg) {
  using F = GF<3, 1>;
  EXPECT_EQ(F::Mul(F{2}, F{2}), F{1});
  EXPECT_EQ(F::Mul(F{2}, F{1}), F{2});
  EXPECT_EQ(F::Neg(F{1}), F{2});
  EXPECT_EQ(F::Neg(F{2}), F{1});
  EXPECT_EQ(F::Neg(F{0}), F{0});
}

TEST(GF3Test, Inverse) {
  using F = GF<3, 1>;
  EXPECT_EQ(F::Inverse(F{1}), F{1});
  EXPECT_EQ(F::Inverse(F{2}), F{2}); // 2*2 = 4 = 1 mod 3
  for (uint8_t a = 1; a < 3; ++a) {
    EXPECT_EQ(F::Mul(F{a}, F::Inverse(F{a})), F::One()) << "a=" << int{a};
  }
}

TEST(GF5Test, Inverse) {
  using F = GF<5, 1>;
  EXPECT_EQ(F::Inverse(F{1}), F{1});
  EXPECT_EQ(F::Inverse(F{2}), F{3}); // 2*3 = 6 = 1 mod 5
  EXPECT_EQ(F::Inverse(F{3}), F{2});
  EXPECT_EQ(F::Inverse(F{4}), F{4}); // 4*4 = 16 = 1 mod 5
  for (uint8_t a = 1; a < 5; ++a) {
    EXPECT_EQ(F::Mul(F{a}, F::Inverse(F{a})), F::One()) << "a=" << int{a};
  }
}

TEST(GF7Test, Inverse) {
  using F = GF<7, 1>;
  for (uint8_t a = 1; a < 7; ++a) {
    EXPECT_EQ(F::Mul(F{a}, F::Inverse(F{a})), F::One()) << "a=" << int{a};
  }
}

TEST(GFConstexprTest, AllOpsAreConstexpr) {
  // Compile-time evaluation pins the constexpr-ness of every op.
  static_assert(GF<3, 1>::Add(1, 2) == GF<3, 1>::Zero());
  static_assert(GF<3, 1>::Sub(1, 2) == GF<3, 1>{2});
  static_assert(GF<3, 1>::Mul(2, 2) == GF<3, 1>::One());
  static_assert(GF<3, 1>::Neg(1) == GF<3, 1>{2});
  static_assert(GF<3, 1>::Inverse(2) == GF<3, 1>{2});
  static_assert(GF<5, 1>::Inverse(2) == GF<5, 1>{3});
  SUCCEED();
}

// Lex-smallest monic irreducible of degree M over F_P, returned as the M low
// coefficients (the M-th is implicitly 1). Reference data: binary extensions
// (q ∈ {4, 8, 16}) match the lex-smallest primitive polynomials tabulated at
// https://www.partow.net/programming/polynomials/index.html. M=1 cases are
// trivially y, and the GF(9) entry y^2 + 1 is the lex-smallest irreducible of
// degree 2 over F_3 (verified by hand: y=0,1,2 are all non-roots).
TEST(GroundIrreducibleTest, AllPrimePowersUpTo16) {
  using gf_internal::GroundIrreducible;
  using A1 = std::array<uint8_t, 1>;
  using A2 = std::array<uint8_t, 2>;
  using A3 = std::array<uint8_t, 3>;
  using A4 = std::array<uint8_t, 4>;
  static_assert(GroundIrreducible<2, 1>() == A1{0});       // q=2:  y
  static_assert(GroundIrreducible<3, 1>() == A1{0});       // q=3:  y
  static_assert(GroundIrreducible<2, 2>() == A2{1, 1});    // q=4:  y^2 + y + 1
  static_assert(GroundIrreducible<5, 1>() == A1{0});       // q=5:  y
  static_assert(GroundIrreducible<7, 1>() == A1{0});       // q=7:  y
  static_assert(GroundIrreducible<2, 3>() == A3{1, 1, 0}); // q=8:  y^3 + y + 1
  static_assert(GroundIrreducible<3, 2>() == A2{1, 0});    // q=9:  y^2 + 1
  static_assert(GroundIrreducible<11, 1>() == A1{0});      // q=11: y
  static_assert(GroundIrreducible<13, 1>() == A1{0});      // q=13: y
  static_assert(GroundIrreducible<2, 4>() ==
                A4{1, 1, 0, 0}); // q=16: y^4 + y + 1
  SUCCEED();
}

// --- M ≥ 2 coverage (lookup-table GF) ---------------------------------------
//
// Elements are F_q indices in [0, q). The bijection index ↔ digit-vector
// matches the natural base-P expansion: index k decodes to digits
// (k mod P, (k/P) mod P, …), i.e. polynomial Σ a_i y^i.

// Verify field axioms over F_q by exhaustive check.
template <int P, int M> void CheckFieldAxioms() {
  using F = GF<P, M>;
  constexpr int q = F::kQ;
  // Additive identity & inverse.
  for (int a = 0; a < q; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Add(fa, F::Zero()), fa) << "a=" << a;
    EXPECT_EQ(F::Add(fa, F::Neg(fa)), F::Zero()) << "a=" << a;
    EXPECT_EQ(F::Sub(fa, fa), F::Zero()) << "a=" << a;
  }
  // Commutativity of Add.
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      F fa{static_cast<uint8_t>(a)};
      F fb{static_cast<uint8_t>(b)};
      EXPECT_EQ(F::Add(fa, fb), F::Add(fb, fa));
    }
  }
  // Multiplicative identity, inverses on nonzero, Inverse(0) = 0.
  EXPECT_EQ(F::Inverse(F::Zero()), F::Zero());
  for (int a = 1; a < q; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Mul(fa, F::One()), fa) << "a=" << a;
    EXPECT_EQ(F::Mul(fa, F::Inverse(fa)), F::One()) << "a=" << a;
  }
  // Multiplicative zero.
  for (int a = 0; a < q; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Mul(fa, F::Zero()), F::Zero()) << "a=" << a;
  }
  // Commutativity of Mul.
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      F fa{static_cast<uint8_t>(a)};
      F fb{static_cast<uint8_t>(b)};
      EXPECT_EQ(F::Mul(fa, fb), F::Mul(fb, fa));
    }
  }
  // Distributivity across all triples.
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      for (int c = 0; c < q; ++c) {
        F fa{static_cast<uint8_t>(a)};
        F fb{static_cast<uint8_t>(b)};
        F fc{static_cast<uint8_t>(c)};
        EXPECT_EQ(F::Mul(fa, F::Add(fb, fc)),
                  F::Add(F::Mul(fa, fb), F::Mul(fa, fc)))
            << "a=" << a << " b=" << b << " c=" << c;
      }
    }
  }
}

TEST(GF4Test, FieldAxioms) { CheckFieldAxioms<2, 2>(); }
TEST(GF8Test, FieldAxioms) { CheckFieldAxioms<2, 3>(); }
TEST(GF9Test, FieldAxioms) { CheckFieldAxioms<3, 2>(); }
TEST(GF16Test, FieldAxioms) { CheckFieldAxioms<2, 4>(); }

// Pow: identity at exp = 0; Fermat at exp = q − 1; matches iterated Mul for
// arbitrary small exponents.
template <int P, int M> void CheckPow() {
  using F = GF<P, M>;
  constexpr int q = F::kQ;
  for (int a = 0; a < q; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Pow(fa, 0), F::One()) << "a=" << a;
  }
  for (int a = 1; a < q; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Pow(fa, q - 1), F::One()) << "a=" << a; // Fermat
  }
  // Cross-check against iterated Mul for small exponents.
  for (int a = 0; a < q; ++a) {
    F fa{static_cast<uint8_t>(a)};
    F acc = F::One();
    for (int e = 0; e < std::min(8, q + 2); ++e) {
      EXPECT_EQ(F::Pow(fa, e), acc) << "a=" << a << " e=" << e;
      acc = F::Mul(acc, fa);
    }
  }
}

TEST(GFPowTest, F2) { CheckPow<2, 1>(); }
TEST(GFPowTest, F3) { CheckPow<3, 1>(); }
TEST(GFPowTest, F5) { CheckPow<5, 1>(); }
TEST(GFPowTest, F7) { CheckPow<7, 1>(); }
TEST(GFPowTest, F4) { CheckPow<2, 2>(); }
TEST(GFPowTest, F8) { CheckPow<2, 3>(); }
TEST(GFPowTest, F9) { CheckPow<3, 2>(); }
TEST(GFPowTest, F16) { CheckPow<2, 4>(); }

// Pow(a, p) is the Frobenius. On F_p it is the identity (Fermat); on a proper
// extension it is a non-trivial field automorphism.
TEST(GFPowTest, FrobeniusOnF4IsNonTrivial) {
  using F = GF<2, 2>;
  // F_4 = F_2[y]/(y^2+y+1). φ(y) = y^2 = y+1, φ(y+1) = y. φ²= id.
  EXPECT_EQ(F::Pow(F{0}, 2), F{0});
  EXPECT_EQ(F::Pow(F{1}, 2), F{1});
  EXPECT_EQ(F::Pow(F{2}, 2), F{3}); // y ↦ y+1
  EXPECT_EQ(F::Pow(F{3}, 2), F{2}); // y+1 ↦ y
  for (int a = 0; a < F::kQ; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Pow(F::Pow(fa, 2), 2), fa) << "a=" << a; // φ² = id
  }
}

TEST(GFPowTest, FrobeniusOnF9IsNonTrivial) {
  using F = GF<3, 2>;
  // F_9 = F_3[y]/(y^2+1). φ(a) = a^3. Fixed points are exactly F_3 ⊂ F_9
  // (indices 0, 1, 2). φ² = id everywhere.
  for (int a = 0; a < 3; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_EQ(F::Pow(fa, 3), fa) << "F_3 fixed: a=" << a;
  }
  for (int a = 3; a < F::kQ; ++a) {
    F fa{static_cast<uint8_t>(a)};
    EXPECT_NE(F::Pow(fa, 3), fa) << "F_9 \\ F_3 moved: a=" << a;
    EXPECT_EQ(F::Pow(F::Pow(fa, 3), 3), fa) << "φ²=id: a=" << a;
  }
}

TEST(GF4Test, KnownProducts) {
  // F_4 = F_2[y]/(y^2 + y + 1). Encoding: 0,1,y,y+1 → 0,1,2,3.
  using F = GF<2, 2>;
  EXPECT_EQ(F::Mul(F{2}, F{2}), F{3}); // y · y = y + 1
  EXPECT_EQ(F::Mul(F{2}, F{3}), F{1}); // y · (y+1) = y^2 + y = 1
  EXPECT_EQ(F::Mul(F{3}, F{3}), F{2}); // (y+1)·(y+1) = y^2 + 1 = y
  EXPECT_EQ(F::Inverse(F{2}), F{3});
  EXPECT_EQ(F::Inverse(F{3}), F{2});
  EXPECT_EQ(F::Add(F{2}, F{3}), F{1});
  EXPECT_EQ(F::Neg(F{3}), F{3}); // F_2 char ⇒ Neg is identity
}

TEST(GFTablesAreConstexpr, F4) {
  using F = GF<2, 2>;
  static_assert(F::Add(2, 3) == F{1});
  static_assert(F::Mul(2, 2) == F{3});
  static_assert(F::Inverse(2) == F{3});
  static_assert(F::Sub(1, 3) == F{2});
  static_assert(F::Neg(3) == F{3});
  SUCCEED();
}

// --- Operator overloads -----------------------------------------------------

TEST(GFOperatorsTest, F4ArithmeticOperators) {
  using F = GF<2, 2>;
  F a{2}, b{3};
  EXPECT_EQ(a + b, F{1});
  EXPECT_EQ(a - b, F{1}); // char 2
  EXPECT_EQ(a * b, F{1});
  EXPECT_EQ(-a, F{2}); // identity in char 2
  EXPECT_EQ(a.Inverse(), F{3});

  F c{2};
  c += b;
  EXPECT_EQ(c, F{1});
  c = F{2};
  c -= b;
  EXPECT_EQ(c, F{1});
  c = F{2};
  c *= b;
  EXPECT_EQ(c, F{1});
}

TEST(GFOperatorsTest, ImplicitFromUint8) {
  using F = GF<3, 1>;
  F a = 2;
  F b = 1;
  EXPECT_EQ(a + b, F::Zero());
}

// --- ToString ---------------------------------------------------------------

TEST(GFToStringTest, F2) {
  using F = GF<2, 1>;
  EXPECT_EQ(F{0}.ToString(), "0");
  EXPECT_EQ(F{1}.ToString(), "1");
}

TEST(GFToStringTest, F3) {
  using F = GF<3, 1>;
  EXPECT_EQ(F{0}.ToString(), "0");
  EXPECT_EQ(F{1}.ToString(), "1");
  EXPECT_EQ(F{2}.ToString(), "2");
}

TEST(GFToStringTest, F13) {
  using F = GF<13, 1>;
  EXPECT_EQ(F{0}.ToString(), "0");
  EXPECT_EQ(F{5}.ToString(), "5");
  EXPECT_EQ(F{9}.ToString(), "9");
  EXPECT_EQ(F{10}.ToString(), "10");
  EXPECT_EQ(F{12}.ToString(), "12");
}

TEST(GFToStringTest, F4) {
  // Encoding: 0,1,y,y+1 → 0,1,2,3; digits low-coef first.
  using F = GF<2, 2>;
  EXPECT_EQ(F{0}.ToString(), "00"); // 0
  EXPECT_EQ(F{1}.ToString(), "10"); // 1
  EXPECT_EQ(F{2}.ToString(), "01"); // y
  EXPECT_EQ(F{3}.ToString(), "11"); // y + 1
}

TEST(GFToStringTest, F8) {
  using F = GF<2, 3>;
  EXPECT_EQ(F{0}.ToString(), "000"); // 0
  EXPECT_EQ(F{1}.ToString(), "100"); // 1
  EXPECT_EQ(F{2}.ToString(), "010"); // y
  EXPECT_EQ(F{4}.ToString(), "001"); // y^2
  EXPECT_EQ(F{7}.ToString(), "111"); // y^2 + y + 1
}

TEST(GFToStringTest, F9) {
  using F = GF<3, 2>;
  EXPECT_EQ(F{0}.ToString(), "00"); // 0
  EXPECT_EQ(F{1}.ToString(), "10"); // 1
  EXPECT_EQ(F{2}.ToString(), "20"); // 2
  EXPECT_EQ(F{3}.ToString(), "01"); // y
  EXPECT_EQ(F{5}.ToString(), "21"); // 2 + y
  EXPECT_EQ(F{8}.ToString(), "22"); // 2 + 2y
}

TEST(GFToStringTest, F16) {
  using F = GF<2, 4>;
  EXPECT_EQ(F{0}.ToString(), "0000");
  EXPECT_EQ(F{1}.ToString(), "1000");
  EXPECT_EQ(F{8}.ToString(), "0001"); // y^3
  EXPECT_EQ(F{15}.ToString(), "1111");
}

TEST(GFOperatorsTest, Layout) {
  static_assert(sizeof(GF<2, 1>) == 1);
  static_assert(sizeof(GF<3, 1>) == 1);
  static_assert(sizeof(GF<2, 2>) == 1);
  static_assert(sizeof(GF<3, 2>) == 1);
  static_assert(sizeof(GF<2, 4>) == 1);
  static_assert(std::is_trivially_copyable_v<GF<2, 1>>);
  static_assert(std::is_trivially_copyable_v<GF<3, 2>>);
  SUCCEED();
}

} // namespace
