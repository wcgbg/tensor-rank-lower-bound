#include "problems/matrix/symmetry.h"

#include <algorithm>
#include <array>
#include <random>
#include <set>

#include <gtest/gtest.h>

#include "core/bit_vec.h"
#include "core/constraints.h"
#include "core/gf.h"
#include "core/gf_vec.h"
#include "core/rank_lower_bound_flatten.h" // FlattenTensorA, CyclicTranspose
#include "core/tensor.h"
#include "problems/matrix/tensor.h"

namespace matrix {
namespace {

// Identity acts trivially; every element is a bijection on functionals whose
// ApplyInverse exactly inverts Apply (over 𝔽_2 the Store canonicalisation is a
// no-op, so the stored matrices satisfy M_inv · M = I).
template <int N0, int N1, int N2> void CheckGroupF2() {
  constexpr int NA = N0 * N1;
  SymmetryGroup<2, 1, N0, N1, N2> g;
  using Vec = GFVec<2, 1, NA>;
  using BV = BitVec<NA>;
  auto make = [](int v) { return Vec{static_cast<BV>(v)}; };
  const int domain = 1 << NA;

  for (int v = 0; v < domain; ++v) {
    EXPECT_EQ(g.store.Apply(g.store.Identity(), make(v)), make(v));
    EXPECT_EQ(g.query.Apply(g.query.Identity(), make(v)), make(v));
  }

  // Elem now carries the matrix data, not the index, so go through At(i).
  for (int t = 0; t < g.store.Size(); ++t) {
    const auto store = g.store.At(t);
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.store.Apply(store, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain) << "store t=" << t;
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.store.ApplyInverse(store, g.store.Apply(store, bv)), bv);
    }
  }

  for (int s = 0; s < g.query.Size(); ++s) {
    const auto query = g.query.At(s);
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.query.Apply(query, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain) << "query s=" << s;
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.query.ApplyInverse(query, g.query.Apply(query, bv)), bv);
    }
  }
}

TEST(SymmetryGroupTest, ActionIsValidForSmallFormats) {
  CheckGroupF2<2, 2, 2>();
  CheckGroupF2<2, 2, 3>();
  CheckGroupF2<2, 3, 3>();
  CheckGroupF2<3, 3, 3>();
}

// Store = |GL(N1, 2)|, Query = (cubic ? 2 : 1) · |GL(N0, 2)|.
// |GL(2,2)| = 6, |GL(3,2)| = 168, |GL(4,2)| = 20160.
TEST(SymmetryGroupTest, GroupSizes) {
  {
    SymmetryGroup<2, 1, 2, 2, 2> g; // cubic
    EXPECT_EQ(g.store.Size(), 6);
    EXPECT_EQ(g.query.Size(), 12);
  }
  {
    SymmetryGroup<2, 1, 3, 3, 3> g; // cubic
    EXPECT_EQ(g.store.Size(), 168);
    EXPECT_EQ(g.query.Size(), 336);
  }
  {
    SymmetryGroup<2, 1, 2, 2, 3> g; // non-cube
    EXPECT_EQ(g.store.Size(), 6);   // GL(2,2)
    EXPECT_EQ(g.query.Size(), 6);   // GL(2,2), no transpose
  }
  {
    SymmetryGroup<2, 1, 2, 3, 4> g; // non-cube
    EXPECT_EQ(g.store.Size(), 168); // GL(3,2)
    EXPECT_EQ(g.query.Size(), 6);   // GL(2,2)
  }
  {
    SymmetryGroup<2, 1, 3, 4, 4> g;   // non-cube
    EXPECT_EQ(g.store.Size(), 20160); // GL(4,2)
    EXPECT_EQ(g.query.Size(), 168);   // GL(3,2)
  }
}

// A random subspace of (𝔽_q^NA)*, returned as its column-reversed RREF basis.
template <int P, int M, int NA>
Constraints<P, M, NA> RandomConstraints(std::mt19937_64 &rng) {
  using Vec = GFVec<P, M, NA>;
  using F = GF<P, M>;
  std::uniform_int_distribution<int> dim_dist(0, NA);
  std::uniform_int_distribution<int> digit_dist(0, F::kQ - 1);
  Constraints<P, M, NA> rows;
  const int target_dim = dim_dist(rng);
  for (int i = 0; i < target_dim; ++i) {
    Vec v{};
    for (int j = 0; j < NA; ++j) {
      v.Set(j, F{static_cast<uint8_t>(digit_dist(rng))});
    }
    rows.push_back(v);
  }
  const int rank = GaussJordanRREF<P, M, NA>(&rows);
  rows.erase(rows.begin(), rows.end() - rank);
  return rows;
}

// The three single-mode flattening ranks of a tensor — each an isomorphism
// invariant. Rotates the axes with CyclicTranspose so one FlattenTensorA covers
// all three positions (non-cube safe).
template <int P, int M, int N0, int N1, int N2>
std::array<int, 3>
ThreeRanks(const Tensor<P, M, N0 * N1, N1 * N2, N2 * N0> &t) {
  constexpr int NA = N0 * N1, NB = N1 * N2, NC = N2 * N0;
  const auto t1 = CyclicTranspose<P, M, NA, NB, NC>(t);
  const auto t2 = CyclicTranspose<P, M, NB, NC, NA>(t1);
  return {
      FlattenTensorA<P, M, NA, NB, NC>(t).Rank(),
      FlattenTensorA<P, M, NB, NC, NA>(t1).Rank(),
      FlattenTensorA<P, M, NC, NA, NB>(t2).Rank(),
  };
}

// Orbit invariance: each (query, store) element acts on the dual A-space and
// lifts to a tensor automorphism of T_mat, so constraining T_mat by R and by
// (query·store)·R must yield isomorphic constrained tensors. We check the
// necessary, cheap invariant that the three flattening ranks agree as a
// multiset (sorted): a pure sandwich element fixes each mode, but the transpose
// element is the ⟨n,n,n⟩ symmetry that swaps the B and C modes, so the ranks
// are preserved only up to that permutation. (Brute-force AreTensorsIsomorphic
// is infeasible here: it enumerates GL(NA)×GL(NB) with NA, NB ≥ 4.)
template <int P, int M, int N0, int N1, int N2>
void CheckConstraintSymmetryInvariance(std::mt19937_64 &rng, int trials) {
  constexpr int NA = N0 * N1, NB = N1 * N2, NC = N2 * N0;
  SymmetryGroup<P, M, N0, N1, N2> g;
  using Tn = Tensor<P, M, NA, NB, NC>;
  const Tn tensor = BuildMulTensor<P, M, N0, N1, N2>();

  for (int trial = 0; trial < trials; ++trial) {
    const Constraints<P, M, NA> constraints = RandomConstraints<P, M, NA>(rng);
    const Tn t0 =
        ApplyConstraintsToTensor<P, M, NA, NB, NC>(constraints, tensor);
    std::array<int, 3> r0 = ThreeRanks<P, M, N0, N1, N2>(t0);
    std::sort(r0.begin(), r0.end());

    for (int s = 0; s < g.query.Size(); ++s) {
      const auto query = g.query.At(s);
      for (int t = 0; t < g.store.Size(); ++t) {
        const auto store = g.store.At(t);
        Constraints<P, M, NA> transformed;
        for (const auto &r : constraints) {
          transformed.push_back(g.query.Apply(query, g.store.Apply(store, r)));
        }
        const int transformed_dim = GaussJordanRREF<P, M, NA>(&transformed);
        EXPECT_EQ(transformed_dim, static_cast<int>(constraints.size()));

        const Tn t1 =
            ApplyConstraintsToTensor<P, M, NA, NB, NC>(transformed, tensor);
        std::array<int, 3> r1 = ThreeRanks<P, M, N0, N1, N2>(t1);
        std::sort(r1.begin(), r1.end());
        EXPECT_EQ(r0, r1) << "N0=" << N0 << " N1=" << N1 << " N2=" << N2
                          << " trial=" << trial << " query.l=" << query.l
                          << " query.t=" << query.transpose
                          << " store=" << store;
      }
    }
  }
}

TEST(SymmetryGroupTest, ConstraintCommutesWithSymmetry) {
  std::mt19937_64 rng;
  CheckConstraintSymmetryInvariance<2, 1, 2, 2, 2>(rng, 30);
  CheckConstraintSymmetryInvariance<2, 1, 2, 2, 3>(rng, 30); // non-cube
  CheckConstraintSymmetryInvariance<2, 1, 2, 3, 4>(rng, 10); // non-cube
  CheckConstraintSymmetryInvariance<2, 1, 3, 3, 3>(rng, 2);  // cubic, big group
}

} // namespace
} // namespace matrix
