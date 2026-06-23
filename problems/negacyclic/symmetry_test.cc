#include "problems/negacyclic/symmetry.h"

#include <array>
#include <random>
#include <set>

#include <gtest/gtest.h>

#include "core/bit_vec.h"
#include "core/constraints.h"
#include "core/gf.h"
#include "core/gf_vec.h"
#include "core/rank_lower_bound_flatten.h"
#include "core/tensor.h"
#include "core/tensor_isomorphism.h"
#include "core/tensor_utils.h"
#include "problems/cyclic/symmetry.h"
#include "problems/negacyclic/tensor.h"

namespace negacyclic {
namespace {

template <int P, int M, int N>
bool SameOneDimSubspace(const GFVec<P, M, N> &a, const GFVec<P, M, N> &b) {
  bool a_zero = true;
  bool b_zero = true;
  for (int i = 0; i < N; ++i) {
    if (a[i].value != 0)
      a_zero = false;
    if (b[i].value != 0)
      b_zero = false;
  }
  if (a_zero != b_zero)
    return false;
  if (a_zero)
    return true;
  int piv = 0;
  while (a[piv].value == 0)
    ++piv;
  if (b[piv].value == 0)
    return false;
  using F = GF<P, M>;
  const F lam = F::Mul(b[piv], F::Inverse(a[piv]));
  for (int i = 0; i < N; ++i) {
    if (!(b[i] == F::Mul(lam, a[i])))
      return false;
  }
  return true;
}

template <int N> void CheckGroupF2() {
  SymmetryGroup<2, 1, N> g;
  using Vec = typename SymmetryGroup<2, 1, N>::Vec;
  using BV = BitVec<N>;
  auto make = [](int v) { return Vec{static_cast<BV>(v)}; };
  const int domain = 1 << N;

  for (int v = 0; v < domain; ++v) {
    EXPECT_EQ(g.store.Apply(g.store.Identity(), make(v)), make(v));
    EXPECT_EQ(g.query.Apply(g.query.Identity(), make(v)), make(v));
  }
  for (int t = 0; t < g.store.Size(); ++t) {
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.store.Apply(t, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain) << "N=" << N << " t=" << t;
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.store.ApplyInverse(t, g.store.Apply(t, bv)), bv);
    }
  }
  for (int s = 0; s < g.query.Size(); ++s) {
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.query.Apply(s, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain) << "N=" << N << " s=" << s;
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.query.ApplyInverse(s, g.query.Apply(s, bv)), bv);
    }
  }
}

TEST(SymmetryGroupTest, ActionIsValidForSmallN) {
  CheckGroupF2<1>();
  CheckGroupF2<2>();
  CheckGroupF2<3>();
  CheckGroupF2<4>();
  CheckGroupF2<5>();
  CheckGroupF2<7>();
}

// Over characteristic 2, x^N + 1 = x^N - 1, so the negacyclic group is *exactly*
// the cyclic group. Construct both and compare sizes for q = 2, 4, 8.
template <int P, int M, int N> void CheckMatchesCyclic() {
  SymmetryGroup<P, M, N> neg;
  cyclic::SymmetryGroup<P, M, N> cyc;
  EXPECT_EQ(neg.store.Size(), cyc.store.Size())
      << "P=" << P << " M=" << M << " N=" << N;
  EXPECT_EQ(neg.query.Size(), cyc.query.Size())
      << "P=" << P << " M=" << M << " N=" << N;
}

TEST(SymmetryGroupTest, MatchesCyclicOverChar2) {
  CheckMatchesCyclic<2, 1, 3>();
  CheckMatchesCyclic<2, 1, 4>();
  CheckMatchesCyclic<2, 1, 5>();
  CheckMatchesCyclic<2, 1, 7>();
  CheckMatchesCyclic<2, 2, 2>(); // q=4
  CheckMatchesCyclic<2, 2, 3>(); // q=4
  CheckMatchesCyclic<2, 3, 2>(); // q=8
}

// Ground-truth sizes for odd prime fields (M=1), brute-forced from the unit /
// auto-image enumeration of R = F_p[x]/(x^N + 1).
TEST(SymmetryGroupTest, StoreQuerySizesOddPrime) {
  auto check = [](int store, int query, auto g) {
    EXPECT_EQ(g.store.Size(), store);
    EXPECT_EQ(g.query.Size(), query);
  };
  check(4, 2, SymmetryGroup<3, 1, 2>{});   // x^2+1 irreducible -> F_9
  check(9, 6, SymmetryGroup<3, 1, 3>{});
  check(32, 8, SymmetryGroup<3, 1, 4>{});  // x^4+1 -> F_9 x F_9
  check(80, 4, SymmetryGroup<3, 1, 5>{});
  check(4, 2, SymmetryGroup<5, 1, 2>{});   // x^2+1 splits -> F_5 x F_5
  check(144, 8, SymmetryGroup<5, 1, 4>{});
  check(8, 2, SymmetryGroup<7, 1, 2>{});
  check(36, 6, SymmetryGroup<7, 1, 3>{});
}

TEST(SymmetryGroupTest, FactorizationIsExact) {
  auto order = [](auto g) { return g.query.Size() * g.store.Size(); };
  EXPECT_EQ(order(SymmetryGroup<3, 1, 2>{}), 8);
  EXPECT_EQ(order(SymmetryGroup<3, 1, 3>{}), 54);
  EXPECT_EQ(order(SymmetryGroup<3, 1, 4>{}), 256);
  EXPECT_EQ(order(SymmetryGroup<5, 1, 2>{}), 8);
  EXPECT_EQ(order(SymmetryGroup<7, 1, 3>{}), 216);
}

template <int P, int M, int N> void CheckFrobeniusInverseIsExact() {
  SymmetryGroup<P, M, N> g;
  using Vec = typename SymmetryGroup<P, M, N>::Vec;
  constexpr uint64_t kQN = IntPow(GF<P, M>::kQ, N);
  for (uint64_t code = 0; code < kQN; ++code) {
    const Vec v = DecodeGFVec<P, M, N>(code);
    for (int r = 0; r < g.query.Size(); ++r) {
      EXPECT_EQ(g.query.ApplyInverse(r, g.query.Apply(r, v)), v)
          << "P=" << P << " M=" << M << " N=" << N << " r=" << r;
    }
  }
}

TEST(SymmetryGroupTest, FrobeniusInverseIsExact) {
  CheckFrobeniusInverseIsExact<3, 1, 2>();
  CheckFrobeniusInverseIsExact<2, 2, 2>(); // q=4
  CheckFrobeniusInverseIsExact<2, 3, 2>(); // q=8
  CheckFrobeniusInverseIsExact<3, 2, 2>(); // q=9
}

template <int P, int M, int N> void CheckProjectiveInverseOverFq() {
  SymmetryGroup<P, M, N> g;
  using Vec = typename SymmetryGroup<P, M, N>::Vec;
  constexpr uint64_t kQN = IntPow(GF<P, M>::kQ, N);
  for (uint64_t code = 0; code < kQN; ++code) {
    const Vec v = DecodeGFVec<P, M, N>(code);
    for (int t = 0; t < g.store.Size(); ++t) {
      const Vec w = g.store.ApplyInverse(t, g.store.Apply(t, v));
      EXPECT_TRUE((SameOneDimSubspace<P, M, N>(v, w)))
          << "P=" << P << " M=" << M << " N=" << N << " t=" << t;
    }
  }
}

TEST(SymmetryGroupTest, InverseIsProjectiveIdentityOverFq) {
  CheckProjectiveInverseOverFq<3, 1, 2>(); // q=3
  CheckProjectiveInverseOverFq<2, 2, 2>(); // q=4
  CheckProjectiveInverseOverFq<3, 2, 2>(); // q=9
}

template <int P, int M, int N>
Constraints<P, M, N> RandomConstraints(std::mt19937_64 &rng) {
  using Vec = GFVec<P, M, N>;
  using F = GF<P, M>;
  std::uniform_int_distribution<int> dim_dist(0, N);
  std::uniform_int_distribution<int> digit_dist(0, F::kQ - 1);
  Constraints<P, M, N> rows;
  const int target_dim = dim_dist(rng);
  for (int i = 0; i < target_dim; ++i) {
    Vec v{};
    for (int j = 0; j < N; ++j) {
      v.Set(j, F{static_cast<uint8_t>(digit_dist(rng))});
    }
    rows.push_back(v);
  }
  const int rank = GaussJordanRREF<P, M, N>(&rows);
  rows.erase(rows.begin(), rows.end() - rank);
  CHECK_EQ(rows.size(), rank);
  CHECK((IsLinearIndependentRREF<P, M, N>(rows)));
  return rows;
}

// Soundness gate: each symmetry element lifts to a tensor automorphism of
// T_neg, so constraining by R and by g·R must give isomorphic tensors.
template <int P, int M, int N>
void CheckConstraintSymmetryIsomorphism(std::mt19937_64 &rng, int fast_trials,
                                        int slow_trials) {
  SymmetryGroup<P, M, N> g;
  using T = Tensor<P, M, N, N, N>;
  const T tensor = BuildMulTensor<P, M, N>();

  auto three_ranks = [](const T &t) -> std::array<int, 3> {
    const auto t1 = CyclicTranspose<P, M, N, N, N>(t);
    const auto t2 = CyclicTranspose<P, M, N, N, N>(t1);
    return {
        FlattenTensorA<P, M, N, N, N>(t).Rank(),
        FlattenTensorA<P, M, N, N, N>(t1).Rank(),
        FlattenTensorA<P, M, N, N, N>(t2).Rank(),
    };
  };

  for (int trial = 0; trial < fast_trials; ++trial) {
    const Constraints<P, M, N> constraints = RandomConstraints<P, M, N>(rng);
    const T t0 = ApplyConstraintsToTensor<P, M, N, N, N>(constraints, tensor);
    const std::array<int, 3> r0 = three_ranks(t0);

    for (int s = 0; s < g.query.Size(); ++s) {
      const auto query = g.query.At(s);
      for (int t = 0; t < g.store.Size(); ++t) {
        const auto store = g.store.At(t);
        Constraints<P, M, N> transformed;
        for (const auto &r : constraints) {
          transformed.push_back(g.query.Apply(query, g.store.Apply(store, r)));
        }
        const int transformed_dim = GaussJordanRREF<P, M, N>(&transformed);
        EXPECT_EQ(transformed_dim, static_cast<int>(constraints.size()));

        const T t1 =
            ApplyConstraintsToTensor<P, M, N, N, N>(transformed, tensor);
        EXPECT_EQ(r0, three_ranks(t1))
            << "P=" << P << " M=" << M << " N=" << N << " trial=" << trial
            << " query=" << query << " store=" << store;

        if (trial < slow_trials) {
          EXPECT_TRUE(AreTensorsIsomorphic(t0, t1))
              << "P=" << P << " M=" << M << " N=" << N << " trial=" << trial
              << " query=" << query << " store=" << store
              << "\n  t0 = " << (TensorToSparseString<P, M, N, N, N>(t0))
              << "\n  t1 = " << (TensorToSparseString<P, M, N, N, N>(t1));
        }
      }
    }
  }
}

TEST(SymmetryGroupTest, ConstraintCommutesWithSymmetryUpToIsomorphism) {
  // slow_trials runs brute-force AreTensorsIsomorphic; restrict it to small
  // (q, N). N ≥ 4, or N = 3 with q ≥ 5, use the rank-only necessary condition.
  std::mt19937_64 rng;
  CheckConstraintSymmetryIsomorphism<2, 1, 3>(rng, 60, 40);
  CheckConstraintSymmetryIsomorphism<3, 1, 2>(rng, 60, 40);
  CheckConstraintSymmetryIsomorphism<3, 1, 3>(rng, 20, 5);  // q=3 N=3
  CheckConstraintSymmetryIsomorphism<3, 1, 4>(rng, 15, 0);  // rank-only
  CheckConstraintSymmetryIsomorphism<2, 2, 2>(rng, 50, 20); // q=4
  CheckConstraintSymmetryIsomorphism<5, 1, 2>(rng, 40, 15); // q=5
  CheckConstraintSymmetryIsomorphism<7, 1, 3>(rng, 15, 0);  // rank-only
  CheckConstraintSymmetryIsomorphism<2, 3, 2>(rng, 30, 5);  // q=8
  CheckConstraintSymmetryIsomorphism<3, 2, 2>(rng, 30, 5);  // q=9
}

} // namespace
} // namespace negacyclic
