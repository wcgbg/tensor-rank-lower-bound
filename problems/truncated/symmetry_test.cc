#include "problems/truncated/symmetry.h"

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
#include "problems/truncated/tensor.h"

namespace truncated {
namespace {

// True iff GFVecs a and b span the same 1-dimensional subspace of 𝔽_q^N.
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
    EXPECT_EQ(g.store.Apply(g.store.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
    EXPECT_EQ(g.query.Apply(g.query.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
  }

  for (int t = 0; t < g.store.Size(); ++t) {
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.store.Apply(t, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain)
        << "N=" << N << " t=" << t;
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.store.ApplyInverse(t, g.store.Apply(t, bv)), bv)
          << "N=" << N << " t=" << t << " v=" << v;
    }
  }

  for (int s = 0; s < g.query.Size(); ++s) {
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.query.Apply(s, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain)
        << "N=" << N << " s=" << s;
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.query.ApplyInverse(s, g.query.Apply(s, bv)), bv)
          << "N=" << N << " s=" << s << " v=" << v;
    }
  }
}

// Closure: composite of any two store elements equals some single element.
// Store = R^*/𝔽_q^* is a subgroup, so this must hold.
template <int N> void CheckClosureF2() {
  SymmetryGroup<2, 1, N> g;
  using Vec = typename SymmetryGroup<2, 1, N>::Vec;
  using BV = BitVec<N>;
  auto make = [](int v) { return Vec{static_cast<BV>(v)}; };
  const int domain = 1 << N;

  for (int s = 0; s < g.store.Size(); ++s) {
    for (int t = 0; t < g.store.Size(); ++t) {
      int found = -1;
      for (int u = 0; u < g.store.Size(); ++u) {
        bool all_match = true;
        for (int v = 0; v < domain && all_match; ++v) {
          const Vec bv = make(v);
          if (g.store.Apply(u, bv) != g.store.Apply(s, g.store.Apply(t, bv))) {
            all_match = false;
          }
        }
        if (all_match) {
          found = u;
          break;
        }
      }
      EXPECT_NE(found, -1) << "N=" << N << " s=" << s << " t=" << t;
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

TEST(SymmetryGroupTest, StoreClosure) {
  CheckClosureF2<1>();
  CheckClosureF2<2>();
  CheckClosureF2<3>();
  CheckClosureF2<4>();
}

// |Store| = |R^*|/(q−1) = q^{N−1} for R = 𝔽_q[x]/(xᴺ) (units are {p_0 ≠ 0}).
TEST(SymmetryGroupTest, StoreSizeMatchesUnitsModScalars) {
  EXPECT_EQ((SymmetryGroup<2, 1, 1>{}.store.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 1, 2>{}.store.Size()), 2);
  EXPECT_EQ((SymmetryGroup<2, 1, 3>{}.store.Size()), 4);
  EXPECT_EQ((SymmetryGroup<2, 1, 4>{}.store.Size()), 8);
  EXPECT_EQ((SymmetryGroup<2, 1, 5>{}.store.Size()), 16);
  EXPECT_EQ((SymmetryGroup<2, 1, 7>{}.store.Size()), 64);

  EXPECT_EQ((SymmetryGroup<3, 1, 2>{}.store.Size()), 3);  // q=3
  EXPECT_EQ((SymmetryGroup<3, 1, 3>{}.store.Size()), 9);  // q=3
  EXPECT_EQ((SymmetryGroup<5, 1, 2>{}.store.Size()), 5);  // q=5
  EXPECT_EQ((SymmetryGroup<2, 2, 2>{}.store.Size()), 4);  // q=4
  EXPECT_EQ((SymmetryGroup<2, 2, 3>{}.store.Size()), 16); // q=4
  EXPECT_EQ((SymmetryGroup<2, 3, 2>{}.store.Size()), 8);  // q=8
  EXPECT_EQ((SymmetryGroup<3, 2, 2>{}.store.Size()), 9);  // q=9
}

// |Query| = M·|Aut(R)| = M·(q−1)·q^{N−2} for N ≥ 2; M for N = 1.
TEST(SymmetryGroupTest, QuerySizeMatchesGaloisTimesAut) {
  EXPECT_EQ((SymmetryGroup<2, 1, 1>{}.query.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 1, 2>{}.query.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 1, 3>{}.query.Size()), 2);
  EXPECT_EQ((SymmetryGroup<2, 1, 4>{}.query.Size()), 4);
  EXPECT_EQ((SymmetryGroup<2, 1, 5>{}.query.Size()), 8);
  EXPECT_EQ((SymmetryGroup<3, 1, 2>{}.query.Size()), 2);  // q=3
  EXPECT_EQ((SymmetryGroup<3, 1, 3>{}.query.Size()), 6);  // q=3
  EXPECT_EQ((SymmetryGroup<5, 1, 2>{}.query.Size()), 4);  // q=5
  EXPECT_EQ((SymmetryGroup<2, 2, 2>{}.query.Size()), 6);  // q=4
  EXPECT_EQ((SymmetryGroup<2, 3, 2>{}.query.Size()), 21); // q=8
  EXPECT_EQ((SymmetryGroup<3, 2, 2>{}.query.Size()), 16); // q=9
}

// Exact split: |Query|·|Store| = |G|.
TEST(SymmetryGroupTest, FactorizationIsExact) {
  auto order = [](auto g) { return g.query.Size() * g.store.Size(); };
  EXPECT_EQ(order(SymmetryGroup<2, 1, 3>{}), 8);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 4>{}), 32);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 5>{}), 128);
  EXPECT_EQ(order(SymmetryGroup<3, 1, 2>{}), 6);    // q=3
  EXPECT_EQ(order(SymmetryGroup<3, 1, 3>{}), 54);   // q=3
  EXPECT_EQ(order(SymmetryGroup<5, 1, 2>{}), 20);   // q=5
  EXPECT_EQ(order(SymmetryGroup<2, 2, 2>{}), 24);   // q=4
  EXPECT_EQ(order(SymmetryGroup<2, 2, 3>{}), 384);  // q=4
  EXPECT_EQ(order(SymmetryGroup<2, 3, 2>{}), 168);  // q=8
  EXPECT_EQ(order(SymmetryGroup<3, 2, 2>{}), 144);  // q=9
}

template <int P, int M, int N> void CheckFrobeniusInverseIsExact() {
  SymmetryGroup<P, M, N> g;
  using Vec = typename SymmetryGroup<P, M, N>::Vec;
  constexpr uint64_t kQN = IntPow(GF<P, M>::kQ, N);
  for (uint64_t code = 0; code < kQN; ++code) {
    const Vec v = DecodeGFVec<P, M, N>(code);
    for (int r = 0; r < g.query.Size(); ++r) {
      EXPECT_EQ(g.query.ApplyInverse(r, g.query.Apply(r, v)), v)
          << "P=" << P << " M=" << M << " N=" << N << " r=" << r
          << " code=" << code;
    }
  }
}

TEST(SymmetryGroupTest, FrobeniusInverseIsExact) {
  CheckFrobeniusInverseIsExact<2, 1, 3>();
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
          << "P=" << P << " M=" << M << " N=" << N << " t=" << t
          << " code=" << code;
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

// Orbit invariance: each symmetry element lifts to a tensor automorphism of
// T_trunc, so constraining T_trunc by R and by g·R must give isomorphic
// tensors. We verify (1) the three flattening ranks agree (necessary, any q),
// and (2) brute-force isomorphism for small q. This is the real soundness gate
// for the symmetry group.
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
  // slow_trials runs the (expensive) brute-force AreTensorsIsomorphic; keep it
  // only for small (q, N) where |GL(N, q)| is tiny. Larger cases use the
  // rank-only necessary condition (slow_trials = 0).
  std::mt19937_64 rng;
  CheckConstraintSymmetryIsomorphism<2, 1, 2>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<2, 1, 3>(rng, 60, 40);
  CheckConstraintSymmetryIsomorphism<2, 1, 4>(rng, 30, 0); // rank-only (N=4)
  CheckConstraintSymmetryIsomorphism<2, 2, 2>(rng, 50, 20); // q=4
  CheckConstraintSymmetryIsomorphism<3, 1, 2>(rng, 60, 40); // q=3
  CheckConstraintSymmetryIsomorphism<5, 1, 2>(rng, 40, 15); // q=5
  CheckConstraintSymmetryIsomorphism<7, 1, 2>(rng, 30, 10); // q=7
  CheckConstraintSymmetryIsomorphism<2, 3, 2>(rng, 30, 5); // q=8
  CheckConstraintSymmetryIsomorphism<3, 2, 2>(rng, 30, 5); // q=9
}

} // namespace
} // namespace truncated
