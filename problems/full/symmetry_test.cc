#include "problems/full/symmetry.h"

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
#include "problems/full/tensor.h"

namespace full {
namespace {

// True iff GFVecs a and b span the same 1-dimensional subspace of 𝔽_q^N (i.e.,
// b is some scalar multiple of a, including both being zero).
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

// For N = 1 every group element collapses to the identity (PGL acts trivially
// on a 1-D coefficient space). For N ≥ 2 over 𝔽_2 the store group is
// PGL(2, 𝔽_2) = GL(2, 𝔽_2) ≅ S₃ of order 6 (P=2 has trivial scalar centre,
// so the PGL quotient coincides with GL).
int ExpectedStoreSizeF2(int n) { return n == 1 ? 1 : 6; }

template <int N> void CheckGroup() {
  SymmetryGroup<2, 1, N> g;
  using Vec = typename SymmetryGroup<2, 1, N>::Vec;
  using BV = BitVec<N>;
  auto make = [](int v) { return Vec{static_cast<BV>(v)}; };
  const int domain = 1 << N; // every subset of coordinates is a functional

  EXPECT_EQ(g.query.Size(), 1) << "N=" << N;
  EXPECT_EQ(g.store.Size(), ExpectedStoreSizeF2(N)) << "N=" << N;

  // The identity index acts trivially.
  for (int v = 0; v < domain; ++v) {
    EXPECT_EQ(g.store.Apply(g.store.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
  }

  for (int t = 0; t < g.store.Size(); ++t) {
    // Each element acts as a bijection on functionals.
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.store.Apply(t, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain)
        << "N=" << N << " t=" << t;

    // store.ApplyInverse really inverts the action.
    for (int v = 0; v < domain; ++v) {
      const Vec bv = make(v);
      EXPECT_EQ(g.store.ApplyInverse(t, g.store.Apply(t, bv)), bv)
          << "N=" << N << " t=" << t << " v=" << v;
    }
  }

  // Closure: the composite of any two elements equals some single element. This
  // confirms {store.Apply(t, ·)} is a genuine group (S₃) under composition.
  // Note this is STRONGER than the orbit map needs (core/symmetry.h asks only
  // for soundness + coverage + inverse-closure, not multiplicative closure); we
  // assert it because this store side does happen to be a full subgroup.
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

  // query.Apply is the identity (trivial query side).
  for (int v = 0; v < domain; ++v) {
    EXPECT_EQ(g.query.Apply(g.query.Identity(), make(v)), make(v));
  }
}

TEST(SymmetryGroupTest, ActionIsValidForSmallN) {
  CheckGroup<1>();
  CheckGroup<2>();
  CheckGroup<3>();
  CheckGroup<4>();
  CheckGroup<5>();
  CheckGroup<6>();
  CheckGroup<7>();
  CheckGroup<8>();
}

// The store group is PGL(2, 𝔽_q) acting via Sym^(N−1) on the dual A-space,
// with order q(q²−1) for N ≥ 2 and 1 for N = 1. Independent of N for N ≥ 2
// because the BFS produces the full projective image from the GL generators.
TEST(SymmetryGroupTest, StoreSizeMatchesPGL) {
  EXPECT_EQ((SymmetryGroup<2, 1, 1>{}.store.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 1, 2>{}.store.Size()), 6);
  EXPECT_EQ((SymmetryGroup<2, 1, 8>{}.store.Size()), 6);

  EXPECT_EQ((SymmetryGroup<3, 1, 1>{}.store.Size()), 1);
  EXPECT_EQ((SymmetryGroup<3, 1, 2>{}.store.Size()), 24);
  EXPECT_EQ((SymmetryGroup<3, 1, 3>{}.store.Size()), 24);
  EXPECT_EQ((SymmetryGroup<3, 1, 5>{}.store.Size()), 24);

  EXPECT_EQ((SymmetryGroup<5, 1, 2>{}.store.Size()), 120);
  EXPECT_EQ((SymmetryGroup<5, 1, 4>{}.store.Size()), 120);

  EXPECT_EQ((SymmetryGroup<7, 1, 2>{}.store.Size()), 336);

  // M ≥ 2: q = P^M.
  EXPECT_EQ((SymmetryGroup<2, 2, 1>{}.store.Size()), 1);   // q=4, N=1
  EXPECT_EQ((SymmetryGroup<2, 2, 2>{}.store.Size()), 60);  // q=4
  EXPECT_EQ((SymmetryGroup<2, 2, 4>{}.store.Size()), 60);  // q=4
  EXPECT_EQ((SymmetryGroup<2, 3, 2>{}.store.Size()), 504); // q=8
  EXPECT_EQ((SymmetryGroup<3, 2, 2>{}.store.Size()), 720); // q=9
  // q = 16: the BFS + inverse loop is O(|elems|² · MatMul) and takes seconds
  // to tens of seconds at |elems| = 4080 here. Kept as a one-off coverage
  // case; deferred perf optimisation tracked separately.
  EXPECT_EQ((SymmetryGroup<2, 4, 2>{}.store.Size()), 4080); // q=16
}

// The query group is Gal(𝔽_q/𝔽_p), cyclic of order M (the Frobenius
// φ: a ↦ a^p). Independent of N. Trivial (size 1) iff M = 1.
TEST(SymmetryGroupTest, QuerySizeMatchesGalois) {
  EXPECT_EQ((SymmetryGroup<2, 1, 2>{}.query.Size()), 1);
  EXPECT_EQ((SymmetryGroup<3, 1, 2>{}.query.Size()), 1);
  EXPECT_EQ((SymmetryGroup<5, 1, 4>{}.query.Size()), 1);
  EXPECT_EQ((SymmetryGroup<2, 2, 2>{}.query.Size()), 2); // q=4
  EXPECT_EQ((SymmetryGroup<2, 3, 2>{}.query.Size()), 3); // q=8
  EXPECT_EQ((SymmetryGroup<3, 2, 2>{}.query.Size()), 2); // q=9
  EXPECT_EQ((SymmetryGroup<2, 4, 2>{}.query.Size()), 4); // q=16
}

// Frobenius φ has order M, so applying φ^r then φ^{M−r} returns v exactly
// (unlike PGL Store where the analogue only holds up to a 1-D subspace).
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
  CheckFrobeniusInverseIsExact<2, 1, 3>(); // trivial M = 1 still exercised
  CheckFrobeniusInverseIsExact<2, 2, 2>(); // q = 4
  CheckFrobeniusInverseIsExact<2, 3, 2>(); // q = 8
  CheckFrobeniusInverseIsExact<3, 2, 2>(); // q = 9
}

// For q > 2 the PGL composition law holds only up to a global scalar on
// individual Vecs: ApplyInverse(t, Apply(t, v)) returns some λ·v with
// λ ∈ 𝔽_q^*. Both vectors span the same 1-D subspace; verify this on every
// vector for a few small (P, M, N).
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
  CheckProjectiveInverseOverFq<3, 1, 2>(); // q = 3
  CheckProjectiveInverseOverFq<5, 1, 2>(); // q = 5
  CheckProjectiveInverseOverFq<2, 2, 2>(); // q = 4
  CheckProjectiveInverseOverFq<3, 2, 2>(); // q = 9
}

// Pin down the concrete generator actions on N = 3, and confirm the
// substitution composed with the reverse has order 3 (the order-3 rotation of
// S₃).
TEST(SymmetryGroupTest, ConcreteActionN3) {
  SymmetryGroup<2, 1, 3> g;
  using Vec = typename SymmetryGroup<2, 1, 3>::Vec;
  auto v3 = [](uint8_t bits) { return Vec{static_cast<BitVec<3>>(bits)}; };

  // Substitution dual G_S: (G_S v)_i = XOR over submasks j ⊆ i of v_j. So the
  // single-coordinate images are: e0 ↦ e0+e1+e2 (every i has 0 as a submask),
  // e1 ↦ e1, e2 ↦ e2 (only i with 1 ⊆ i is i=1; only i with 2 ⊆ i is i=2).
  auto is_substitution = [&](int t) {
    return g.store.Apply(t, v3(0b001)) == v3(0b111) &&
           g.store.Apply(t, v3(0b010)) == v3(0b010) &&
           g.store.Apply(t, v3(0b100)) == v3(0b100);
  };
  // Reverse dual G_σ: bit i ↦ bit N−1−i, i.e. e0 ↔ e2, e1 fixed.
  auto is_reverse = [&](int t) {
    return g.store.Apply(t, v3(0b001)) == v3(0b100) &&
           g.store.Apply(t, v3(0b010)) == v3(0b010) &&
           g.store.Apply(t, v3(0b100)) == v3(0b001);
  };

  int sub = -1;
  int rev = -1;
  for (int t = 0; t < g.store.Size(); ++t) {
    if (is_substitution(t))
      sub = t;
    if (is_reverse(t))
      rev = t;
  }
  ASSERT_NE(sub, -1);
  ASSERT_NE(rev, -1);

  // ρ = S∘σ has order 3: ρ³ = identity, ρ ≠ id, ρ² ≠ id, acting on functionals.
  const int domain = 1 << 3;
  auto apply_rho = [&](Vec v) {
    return g.store.Apply(sub, g.store.Apply(rev, v));
  };
  bool rho_is_id = true;
  bool rho2_is_id = true;
  for (int v = 0; v < domain; ++v) {
    const Vec bv = v3(static_cast<uint8_t>(v));
    const Vec r1 = apply_rho(bv);
    const Vec r2 = apply_rho(r1);
    const Vec r3 = apply_rho(r2);
    if (!(r1 == bv))
      rho_is_id = false;
    if (!(r2 == bv))
      rho2_is_id = false;
    EXPECT_EQ(r3, bv) << "v=" << v; // order divides 3
  }
  EXPECT_FALSE(rho_is_id);
  EXPECT_FALSE(rho2_is_id);
}

// A random subspace of (𝔽_q^N)*, returned as its column-reversed RREF basis
// (the canonical form GaussJordanRREF produces). The dimension is random,
// each coordinate is a uniform digit in [0, q) where q = P^M.
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

// Orbit invariance: a symmetry element g acts on the dual A-space (store.Apply)
// and lifts to a tensor automorphism of T_full. Hence constraining T_full by a
// subspace R and by the transported subspace g·R must give isomorphic tensors.
//
// We verify two invariants:
//   (1) the three flattening ranks (one per mode) agree. This is a *necessary*
//       condition for tensor isomorphism — every invertible (A, B, C)
//       preserves each flattening's row span and hence its rank — and is
//       cheap (one Gaussian elimination per mode). Works for any (P, M).
//   (2) for q ∈ {2, 3} (i.e. (P ≤ 3, M = 1)), the brute-force
//       AreTensorsIsomorphic query. This is O(|GL(NA, 𝔽_q)|·|GL(NB,
//       𝔽_q)|·|GL(NC, 𝔽_q)|), which already costs ~10⁷ tensor multiplications
//       at (P=3, N=2) and ~10¹¹ at (P=5, N=2); only the smallest q fits in a
//       unit-test budget.
template <int P, int M, int N>
void CheckConstraintSymmetryIsomorphism(std::mt19937_64 &rng, int fast_trials,
                                        int slow_trials) {
  SymmetryGroup<P, M, N> g;
  using T = Tensor<P, M, N, N, 2 * N - 1>;
  const T tensor = BuildMulTensor<P, M, N>();

  // Flatten ranks in each of the three modes — invariants of tensor
  // isomorphism. Uses CyclicTranspose to bring each axis to position A.
  auto three_ranks = [](const T &t) -> std::array<int, 3> {
    const auto t1 = CyclicTranspose<P, M, N, N, 2 * N - 1>(t);
    const auto t2 = CyclicTranspose<P, M, N, 2 * N - 1, N>(t1);
    return {
        FlattenTensorA<P, M, N, N, 2 * N - 1>(t).Rank(),
        FlattenTensorA<P, M, N, 2 * N - 1, N>(t1).Rank(),
        FlattenTensorA<P, M, 2 * N - 1, N, N>(t2).Rank(),
    };
  };

  for (int trial = 0; trial < fast_trials; ++trial) {
    const Constraints<P, M, N> constraints = RandomConstraints<P, M, N>(rng);
    const T t0 =
        ApplyConstraintsToTensor<P, M, N, N, 2 * N - 1>(constraints, tensor);
    const std::array<int, 3> r0 = three_ranks(t0);

    // Iterate the full Query × Store product. Query applied after Store,
    // matching the framework's convention (core/orbit_enumerator_slow.h:
    // Transform). For M = 1 the query loop is a single identity iteration
    // and this collapses to the original Store-only sweep.
    for (int s = 0; s < g.query.Size(); ++s) {
      const auto query = g.query.At(s);
      for (int t = 0; t < g.store.Size(); ++t) {
        const auto store = g.store.At(t);

        // Transport the subspace by g, then renormalise to RREF.
        Constraints<P, M, N> transformed;
        for (const auto &r : constraints) {
          transformed.push_back(g.query.Apply(query, g.store.Apply(store, r)));
        }
        const int transformed_dim = GaussJordanRREF<P, M, N>(&transformed);
        EXPECT_EQ(transformed_dim, static_cast<int>(constraints.size()));

        const T t1 = ApplyConstraintsToTensor<P, M, N, N, 2 * N - 1>(
            transformed, tensor);
        const std::array<int, 3> r1 = three_ranks(t1);
        EXPECT_EQ(r0, r1) << "P=" << P << " M=" << M << " N=" << N
                          << " trial=" << trial << " query=" << query
                          << " store=" << store;

        if (trial < slow_trials) {
          EXPECT_TRUE(AreTensorsIsomorphic(t0, t1))
              << "P=" << P << " M=" << M << " N=" << N << " trial=" << trial
              << " query=" << query << " store=" << store << "\n  t0 = "
              << (TensorToSparseString<P, M, N, N, 2 * N - 1>(t0))
              << "\n  t1 = "
              << (TensorToSparseString<P, M, N, N, 2 * N - 1>(t1));
        }
      }
    }
  }
}

TEST(SymmetryGroupTest, ConstraintCommutesWithSymmetryUpToIsomorphism) {
  std::mt19937_64 rng;
  CheckConstraintSymmetryIsomorphism<2, 1, 2>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<3, 1, 2>(rng, 100, 100);
  CheckConstraintSymmetryIsomorphism<2, 2, 2>(rng, 100, 10); // q=4
  CheckConstraintSymmetryIsomorphism<5, 1, 2>(rng, 100, 10);
  CheckConstraintSymmetryIsomorphism<7, 1, 2>(rng, 100, 2);
  CheckConstraintSymmetryIsomorphism<2, 3, 2>(rng, 100, 2); // q=8
  CheckConstraintSymmetryIsomorphism<3, 2, 2>(rng, 100, 1); // q=9
}

} // namespace
} // namespace full
