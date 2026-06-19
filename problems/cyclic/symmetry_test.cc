#include "problems/cyclic/symmetry.h"

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
#include "problems/cyclic/tensor.h"

namespace cyclic {
namespace {

// True iff GFVecs a and b span the same 1-dimensional subspace of 𝔽_q^N (b is
// some scalar multiple of a, both zero counted as the same 0-D subspace).
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

  // Both sides' identity acts trivially.
  for (int v = 0; v < domain; ++v) {
    EXPECT_EQ(g.store.Apply(g.store.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
    EXPECT_EQ(g.query.Apply(g.query.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
  }

  // Each store element acts as a bijection on functionals.
  for (int t = 0; t < g.store.Size(); ++t) {
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.store.Apply(t, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain)
        << "N=" << N << " t=" << t;

    // store.ApplyInverse exactly inverts the action over 𝔽_2 (canonicalisation
    // is a no-op, so the stored matrices satisfy M_inv · M = I).
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.store.ApplyInverse(t, g.store.Apply(t, bv)), bv)
          << "N=" << N << " t=" << t << " v=" << v;
    }
  }

  // Each query element (a ring auto, for M = 1) is likewise a bijection whose
  // ApplyInverse exactly inverts Apply (literal inverse matrices).
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

// Closure: the composite of any two store elements equals some single element.
// O(|store|³ · 2^N), only feasible for small N where |store| stays tiny.
// Stronger than the orbit map needs (core/symmetry.h requires only soundness +
// coverage + inverse-closure); asserted because this store side is a subgroup.
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
}

// Store now counts R^*/𝔽_q^* (the multiplications). Closed-form sizes
// |R^*|/(q−1); see the header comment for the per-N derivations.
TEST(SymmetryGroupTest, StoreSizeMatchesUnitsModScalars) {
  // (P=2, M=1) — q = 2, so |Store| = |R^*|.
  EXPECT_EQ((SymmetryGroup<2, 1, 1>{}.store.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 1, 2>{}.store.Size()), 2);
  EXPECT_EQ((SymmetryGroup<2, 1, 3>{}.store.Size()), 3);
  EXPECT_EQ((SymmetryGroup<2, 1, 4>{}.store.Size()), 8);
  EXPECT_EQ((SymmetryGroup<2, 1, 5>{}.store.Size()), 15);
  EXPECT_EQ((SymmetryGroup<2, 1, 7>{}.store.Size()), 49);
  EXPECT_EQ((SymmetryGroup<2, 1, 9>{}.store.Size()), 189);
  EXPECT_EQ((SymmetryGroup<2, 1, 11>{}.store.Size()), 1023);
  EXPECT_EQ((SymmetryGroup<2, 1, 13>{}.store.Size()), 4095);

  // q > 2.
  EXPECT_EQ((SymmetryGroup<2, 2, 1>{}.store.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 2, 2>{}.store.Size()), 4);
  EXPECT_EQ((SymmetryGroup<2, 2, 3>{}.store.Size()), 9);
  EXPECT_EQ((SymmetryGroup<3, 1, 1>{}.store.Size()), 1);
  EXPECT_EQ((SymmetryGroup<3, 1, 2>{}.store.Size()), 2);
}

// Query now counts Gal × Aut(R), size M·|Aut(R)|.
TEST(SymmetryGroupTest, QuerySizeMatchesGaloisTimesAut) {
  EXPECT_EQ((SymmetryGroup<2, 1, 3>{}.query.Size()), 2);
  EXPECT_EQ((SymmetryGroup<3, 1, 2>{}.query.Size()), 2);
  EXPECT_EQ((SymmetryGroup<5, 1, 2>{}.query.Size()), 2);
  EXPECT_EQ((SymmetryGroup<2, 2, 2>{}.query.Size()), 6);  // q=4
  EXPECT_EQ((SymmetryGroup<2, 3, 2>{}.query.Size()), 21); // q=8
  EXPECT_EQ((SymmetryGroup<3, 2, 2>{}.query.Size()), 4);  // q=9
}

// The split is exact: |Query|·|Store| = |G|, the order it had when the entire
// (R^* ⋊ Aut(R))/𝔽_q^* factor lived on Store and Gal lived on Query.
TEST(SymmetryGroupTest, FactorizationIsExact) {
  auto order = [](auto g) { return g.query.Size() * g.store.Size(); };
  EXPECT_EQ(order(SymmetryGroup<2, 1, 3>{}), 6);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 4>{}), 32);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 5>{}), 60);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 7>{}), 882);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 9>{}), 2268);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 11>{}), 10230);
  EXPECT_EQ(order(SymmetryGroup<2, 1, 13>{}), 49140);
  EXPECT_EQ(order(SymmetryGroup<2, 2, 2>{}), 24);  // q=4
  EXPECT_EQ(order(SymmetryGroup<2, 2, 3>{}), 108); // q=4
  EXPECT_EQ(order(SymmetryGroup<2, 3, 2>{}), 168); // q=8
  EXPECT_EQ(order(SymmetryGroup<3, 1, 2>{}), 4);   // q=3
  EXPECT_EQ(order(SymmetryGroup<3, 2, 2>{}), 32);  // q=9
  EXPECT_EQ(order(SymmetryGroup<5, 1, 2>{}), 8);   // q=5
}

// Every query element (auto a, Galois power r) inverts exactly: the literal
// dual matrices give M_a⁻¹·M_a = I and φ^r ∘ φ^(M−r) = id, so
// ApplyInverse(e, Apply(e, v)) = v on the nose for every vector.
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
  CheckFrobeniusInverseIsExact<2, 1, 3>(); // M=1 trivial branch
  CheckFrobeniusInverseIsExact<2, 2, 2>(); // q=4
  CheckFrobeniusInverseIsExact<2, 3, 2>(); // q=8
  CheckFrobeniusInverseIsExact<3, 2, 2>(); // q=9
}

// For q > 2 the store composition holds only up to a global 𝔽_q^* scalar on
// individual Vecs: ApplyInverse(t, Apply(t, v)) = λ·v with λ ∈ 𝔽_q^*. Both
// vectors span the same 1-D subspace.
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

// A random subspace of (𝔽_q^N)*, returned as its column-reversed RREF basis.
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

// Orbit invariance: each symmetry element g acts on the dual A-space and lifts
// to a tensor automorphism of T_cyc, so constraining T_cyc by R and by g·R
// must give isomorphic tensors. We verify two invariants:
//   (1) the three flattening ranks (one per mode) agree — necessary, cheap,
//       works for any (P, M).
//   (2) brute-force AreTensorsIsomorphic query for small q (only feasible
//       at q = 2 here, with N ≤ 3).
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
  std::mt19937_64 rng;
  CheckConstraintSymmetryIsomorphism<2, 1, 2>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<2, 1, 3>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<2, 2, 2>(rng, 100, 100); // q=4: rank-only
  CheckConstraintSymmetryIsomorphism<5, 1, 2>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<7, 1, 2>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<2, 3, 2>(rng, 100, 10); // q=8
  CheckConstraintSymmetryIsomorphism<3, 2, 2>(rng, 100, 10); // q=9
}

} // namespace
} // namespace cyclic
