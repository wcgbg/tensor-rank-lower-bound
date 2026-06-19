#include "problems/extfield/symmetry.h"

#include <array>
#include <random>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "core/bit_vec.h"
#include "core/constraints.h"
#include "core/tensor.h"
#include "core/tensor_isomorphism.h"
#include "core/tensor_utils.h"
#include "problems/extfield/irreducibles.h"
#include "problems/extfield/tensor.h"

namespace extfield {
namespace {

template <int N> BitVec<N> Modulus() {
  return static_cast<BitVec<N>>(IrreduciblePolyBits<N>().to_ullong());
}

// FieldMul against hand-computed products.
TEST(FieldMulTest, KnownProducts) {
  // GF(8) = 𝔽₂[x]/(x³+x+1): x · x² = x³ = x+1.
  const uint8_t p3 = Modulus<3>();
  EXPECT_EQ(FieldMul<3>(uint8_t{0b010}, uint8_t{0b100}, p3), uint8_t{0b011});
  // x² · x² = x⁴ = x²+x.
  EXPECT_EQ(FieldMul<3>(uint8_t{0b100}, uint8_t{0b100}, p3), uint8_t{0b110});
  // 1 is the multiplicative identity.
  EXPECT_EQ(FieldMul<3>(uint8_t{0b101}, uint8_t{0b001}, p3), uint8_t{0b101});

  // GF(256), AES field 𝔽₂[x]/(x⁸+x⁴+x³+x+1): the worked Wikipedia example
  // (x⁶+x⁴+x+1)·(x⁷+x⁶+x³+x) = 1.
  const uint8_t p8 = Modulus<8>();
  EXPECT_EQ(FieldMul<8>(uint8_t{0b01010011}, uint8_t{0b11001010}, p8),
            uint8_t{0b00000001});
}

// FieldInverse(a) · a == 1 for every nonzero a, at a few small N.
template <int N> void CheckFieldInverse() {
  using BV = BitVec<N>;
  const BV p = Modulus<N>();
  for (int a = 1; a < (1 << N); ++a) {
    const BV av = static_cast<BV>(a);
    const BV inv = FieldInverse<N>(av, p);
    EXPECT_EQ(FieldMul<N>(av, inv, p), BV{1}) << "N=" << N << " a=" << a;
  }
}

TEST(FieldInverseTest, InvertsEveryNonzeroElement) {
  CheckFieldInverse<1>();
  CheckFieldInverse<2>();
  CheckFieldInverse<3>();
  CheckFieldInverse<4>();
  CheckFieldInverse<5>();
  CheckFieldInverse<6>();
  CheckFieldInverse<7>();
  CheckFieldInverse<8>();
}

template <int N> void CheckGroup() {
  SymmetryGroup<2, 1, N> g;
  using Vec = typename SymmetryGroup<2, 1, N>::Vec;
  using BV = BitVec<N>;
  auto make = [](int v) { return Vec{static_cast<BV>(v)}; };
  const int domain = 1 << N; // every subset of coordinates is a functional

  EXPECT_EQ(g.query.Size(), N) << "N=" << N;
  EXPECT_EQ(g.store.Size(), (1 << N) - 1) << "N=" << N;

  // The identity indices act trivially.
  for (int v = 0; v < domain; ++v) {
    EXPECT_EQ(g.query.Apply(g.query.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
    EXPECT_EQ(g.store.Apply(g.store.Identity(), make(v)), make(v))
        << "N=" << N << " v=" << v;
  }

  for (int t = 0; t < g.store.Size(); ++t) {
    // Each store element (scaling) acts as a bijection on functionals.
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

  // Each query (Frobenius power) element is a bijection too.
  for (int s = 0; s < g.query.Size(); ++s) {
    std::set<BV> image;
    for (int v = 0; v < domain; ++v) {
      image.insert(g.query.Apply(s, make(v)).data);
    }
    EXPECT_EQ(static_cast<int>(image.size()), domain)
        << "N=" << N << " s=" << s;
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

// Every group element is uniquely (query⁻¹ · store). We identify a composite
// by its action on the N basis functionals and check there are exactly N·(2ᴺ−1)
// distinct maps — the full group order, with no collisions. Uniqueness is more
// than the orbit map needs (core/symmetry.h requires only that Query⁻¹ · Store
// *cover* the group, not that the factorisation be unique); we assert it here
// because this implementation's exact subgroup factorisation provides it, which
// keeps |G| = query.Size()·store.Size() and one witness per orbit.
template <int N> void CheckFactorizationIsUnique() {
  SymmetryGroup<2, 1, N> g;
  using Vec = typename SymmetryGroup<2, 1, N>::Vec;
  using BV = BitVec<N>;

  std::set<std::array<BV, N>> distinct;
  for (int s = 0; s < g.query.Size(); ++s) {
    for (int t = 0; t < g.store.Size(); ++t) {
      std::array<BV, N> signature{};
      for (int i = 0; i < N; ++i) {
        const Vec ei{static_cast<BV>(BV{1} << i)};
        // Framework order: store first, then query — Apply(g, v) =
        // query.Apply(s, store.Apply(t, v)), matching
        // OrbitEnumeratorSlow::Transform.
        signature[i] = g.query.Apply(s, g.store.Apply(t, ei)).data;
      }
      distinct.insert(signature);
    }
  }
  EXPECT_EQ(static_cast<int>(distinct.size()), N * ((1 << N) - 1)) << "N=" << N;
}

TEST(SymmetryGroupTest, QueryStoreFactorizationIsUnique) {
  CheckFactorizationIsUnique<1>();
  CheckFactorizationIsUnique<2>();
  CheckFactorizationIsUnique<3>();
  CheckFactorizationIsUnique<4>();
  CheckFactorizationIsUnique<5>();
  CheckFactorizationIsUnique<6>();
  CheckFactorizationIsUnique<7>();
  CheckFactorizationIsUnique<8>();
}

// Pin down the concrete action on GF(8) = 𝔽₂[x]/(x³+x+1).
TEST(SymmetryGroupTest, ConcreteActionN3) {
  SymmetryGroup<2, 1, 3> g;
  using Vec = typename SymmetryGroup<2, 1, 3>::Vec;
  using BV3 = BitVec<3>;
  auto make = [](uint8_t bits) { return Vec{static_cast<BV3>(bits)}; };
  const uint8_t domain = 1 << 3;

  // Frobenius (query index 1) has order 3 on the dual space.
  auto apply_phi = [&](Vec v) { return g.query.Apply(1, v); };
  bool phi_is_id = true;
  bool phi2_is_id = true;
  for (int v = 0; v < domain; ++v) {
    const Vec bv = make(static_cast<uint8_t>(v));
    const Vec r1 = apply_phi(bv);
    const Vec r2 = apply_phi(r1);
    const Vec r3 = apply_phi(r2);
    if (!(r1 == bv))
      phi_is_id = false;
    if (!(r2 == bv))
      phi2_is_id = false;
    EXPECT_EQ(r3, bv) << "v=" << v; // order divides 3
  }
  EXPECT_FALSE(phi_is_id);
  EXPECT_FALSE(phi2_is_id);

  // The scaling subgroup is cyclic of order 7. Store index i denotes scaling
  // by the field element α = i+1, and scaling_α ∘ scaling_β = scaling_{αβ}, so
  // the powers of α = x (index 1) are scalings by x^p. Their distinct dual
  // actions must cover all 7 group elements (x is a generator of GF(8)^*).
  const uint8_t p3 = Modulus<3>();
  std::set<std::array<BV3, 3>> powers;
  uint8_t alpha = 1; // x^0
  for (int p = 0; p < g.store.Size(); ++p) {
    const int idx = static_cast<int>(alpha) - 1; // scaling by α
    std::array<BV3, 3> sig{};
    for (int i = 0; i < 3; ++i) {
      sig[i] = g.store.Apply(idx, make(static_cast<uint8_t>(1 << i))).data;
    }
    powers.insert(sig);
    alpha = FieldMul<3>(alpha, uint8_t{0b010}, p3); // α ← α·x
  }
  EXPECT_EQ(static_cast<int>(powers.size()), g.store.Size());
}

// A random subspace of (𝔽₂ᴺ)*, returned as its column-reversed RREF basis. The
// dimension is itself random in [0, N].
template <int N> Constraints<2, 1, N> RandomConstraints(std::mt19937_64 &rng) {
  using Vec = GFVec<2, 1, N>;
  using BV = BitVec<N>;
  std::uniform_int_distribution<int> dim_dist(0, N);
  std::uniform_int_distribution<int> vec_dist(1, (1 << N) - 1);
  Constraints<2, 1, N> rows;
  int dim = dim_dist(rng);
  for (int i = 0; i < dim; ++i) {
    rows.push_back(Vec{static_cast<BV>(vec_dist(rng))});
  }
  dim = GaussJordanRREF<2, 1, N>(&rows);
  rows.erase(rows.begin(), rows.end() - dim);
  CHECK_EQ(rows.size(), dim);
  CHECK((IsLinearIndependentRREF<2, 1, N>(rows)));
  return rows;
}

// Orbit invariance: a symmetry element acts on the dual A-space (query/store)
// and lifts to a tensor automorphism of T_ext. Constraining T_ext by a subspace
// R and by its transported image must give isomorphic tensors. (Brute-force
// isomorphism is GL(N)³, so keep N tiny.)
template <int N>
void CheckConstraintSymmetryIsomorphism(std::mt19937_64 &rng, int trials) {
  SymmetryGroup<2, 1, N> g;
  using T = Tensor<2, 1, N, N, N>;
  const T tensor = BuildMulTensor<N>(IrreduciblePolyBits<N>());
  std::uniform_int_distribution<int> search_dist(0, g.query.Size() - 1);
  std::uniform_int_distribution<int> store_dist(0, g.store.Size() - 1);

  for (int trial = 0; trial < trials; ++trial) {
    const Constraints<2, 1, N> constraints = RandomConstraints<N>(rng);
    const auto query = g.query.At(search_dist(rng));
    const auto store = g.store.At(store_dist(rng));

    // Transport the subspace by query then store, then renormalise to RREF.
    Constraints<2, 1, N> transformed;
    for (const auto &r : constraints) {
      transformed.push_back(g.store.Apply(store, g.query.Apply(query, r)));
    }
    int transformed_dim = GaussJordanRREF<2, 1, N>(&transformed);
    EXPECT_EQ(transformed_dim, static_cast<int>(constraints.size()));

    const T t0 = ApplyConstraintsToTensor<2, 1, N, N, N>(constraints, tensor);
    const T t1 = ApplyConstraintsToTensor<2, 1, N, N, N>(transformed, tensor);

    EXPECT_TRUE(AreTensorsIsomorphic(t0, t1))
        << "N=" << N << " trial=" << trial << " query=" << query
        << " store=" << store
        << "\n  t0 = " << (TensorToSparseString<2, 1, N, N, N>(t0))
        << "\n  t1 = " << (TensorToSparseString<2, 1, N, N, N>(t1));
  }
}

TEST(SymmetryGroupTest, ConstraintCommutesWithSymmetryUpToIsomorphism) {
  std::mt19937_64 rng;
  CheckConstraintSymmetryIsomorphism<2>(rng, 40);
  CheckConstraintSymmetryIsomorphism<3>(rng, 20);
}

// --- General (P, M, N) coverage, including the M ≥ 2 base field 𝔽_q ---

// True iff u and v span the same 1-D subspace (u = c·v for a scalar c), or both
// are zero. The store-side inverse is only a *projective* inverse for q > 2
// (the scalar centre 𝔽_q^* was quotiented out), so its round trip lands on the
// same line, not the same vector.
template <int P, int M, int N>
bool SameLine(GFVec<P, M, N> u, GFVec<P, M, N> v) {
  using Ops = GF<P, M>;
  const bool uz = u.IsZero();
  const bool vz = v.IsZero();
  if (uz || vz) {
    return uz && vz;
  }
  const int idx = v.LeadingNonzeroIdx();
  const Ops c = Ops::Mul(u[idx], Ops::Inverse(v[idx]));
  return u == (c * v);
}

template <int P, int M, int N> void CheckGroupGeneric() {
  SymmetryGroup<P, M, N> g;
  using Vec = GFVec<P, M, N>;
  constexpr int kQ = GF<P, M>::kQ;
  const uint64_t qn = IntPow(kQ, N);
  const int store_expected = static_cast<int>((qn - 1) / (kQ - 1));

  EXPECT_EQ(g.query.Size(), M * N) << "P=" << P << " M=" << M << " N=" << N;
  EXPECT_EQ(g.store.Size(), store_expected)
      << "P=" << P << " M=" << M << " N=" << N;

  // The identity indices act trivially.
  for (uint64_t i = 0; i < qn; ++i) {
    const Vec v = DecodeGFVec<P, M, N>(i);
    EXPECT_EQ(g.query.Apply(g.query.Identity(), v), v);
    EXPECT_EQ(g.store.Apply(g.store.Identity(), v), v);
  }

  // Store: each scaling dual is an invertible linear map (a bijection on every
  // 𝔽_q^N vector); ApplyInverse is a projective inverse (round trip → same
  // line).
  for (int s = 0; s < g.store.Size(); ++s) {
    std::set<uint64_t> image;
    for (uint64_t i = 0; i < qn; ++i) {
      image.insert(
          EncodeGFVec<P, M, N>(g.store.Apply(s, DecodeGFVec<P, M, N>(i))));
    }
    EXPECT_EQ(static_cast<uint64_t>(image.size()), qn) << "store s=" << s;
    for (uint64_t i = 0; i < qn; ++i) {
      const Vec v = DecodeGFVec<P, M, N>(i);
      EXPECT_TRUE(
          (SameLine<P, M, N>(g.store.ApplyInverse(s, g.store.Apply(s, v)), v)))
          << "store s=" << s << " i=" << i;
    }
  }

  // Query: each field automorphism dual is a bijection, with an exact inverse.
  for (int k = 0; k < g.query.Size(); ++k) {
    std::set<uint64_t> image;
    for (uint64_t i = 0; i < qn; ++i) {
      image.insert(
          EncodeGFVec<P, M, N>(g.query.Apply(k, DecodeGFVec<P, M, N>(i))));
    }
    EXPECT_EQ(static_cast<uint64_t>(image.size()), qn) << "query k=" << k;
    for (uint64_t i = 0; i < qn; ++i) {
      const Vec v = DecodeGFVec<P, M, N>(i);
      EXPECT_EQ(g.query.ApplyInverse(k, g.query.Apply(k, v)), v)
          << "query k=" << k << " i=" << i;
    }
  }
}

// Every (query, store) pair gives a distinct composite map (the factorization
// is unique), so the count of distinct maps is query.Size()·store.Size() =
// M·N·(q^N−1)/(q−1). The map is compared by its full action on all q^N
// functionals (semilinear maps are not pinned down by a basis alone), in the
// framework order query.ApplyInverse(s, store.Apply(t, ·)). (Uniqueness is a
// property of this exact subgroup factorisation, not a framework requirement —
// the orbit map needs only that Query⁻¹ · Store covers the group; see
// core/symmetry.h.)
template <int P, int M, int N> void CheckFactorizationGeneric() {
  SymmetryGroup<P, M, N> g;
  const uint64_t qn = IntPow(GF<P, M>::kQ, N);
  std::set<std::vector<uint64_t>> distinct;
  for (int s = 0; s < g.query.Size(); ++s) {
    for (int t = 0; t < g.store.Size(); ++t) {
      std::vector<uint64_t> sig(qn);
      for (uint64_t i = 0; i < qn; ++i) {
        sig[i] = EncodeGFVec<P, M, N>(
            g.query.Apply(s, g.store.Apply(t, DecodeGFVec<P, M, N>(i))));
      }
      distinct.insert(std::move(sig));
    }
  }
  EXPECT_EQ(static_cast<int>(distinct.size()), g.query.Size() * g.store.Size())
      << "P=" << P << " M=" << M << " N=" << N;
}

template <int P, int M, int N>
Constraints<P, M, N> RandomConstraintsGeneric(std::mt19937_64 &rng) {
  const uint64_t qn = IntPow(GF<P, M>::kQ, N);
  std::uniform_int_distribution<int> dim_dist(0, N);
  std::uniform_int_distribution<uint64_t> vec_dist(1, qn - 1);
  Constraints<P, M, N> rows;
  const int dim = dim_dist(rng);
  for (int i = 0; i < dim; ++i) {
    rows.push_back(DecodeGFVec<P, M, N>(vec_dist(rng)));
  }
  const int rank = GaussJordanRREF<P, M, N>(&rows);
  rows.erase(rows.begin(), rows.end() - rank);
  return rows;
}

// A symmetry element transports a constraint subspace (query after store, the
// framework order); the constrained tensor of the image must be isomorphic to
// the original. Brute-force isomorphism is GL(N, q)², so keep this to N = 2.
template <int P, int M, int N>
void CheckConstraintSymmetryIsoGeneric(std::mt19937_64 &rng, int trials) {
  SymmetryGroup<P, M, N> g;
  using T = Tensor<P, M, N, N, N>;
  const T tensor = BuildMulTensor<P, M, N>(IrreduciblePolyCoeffs<P, M, N>());
  std::uniform_int_distribution<int> search_dist(0, g.query.Size() - 1);
  std::uniform_int_distribution<int> store_dist(0, g.store.Size() - 1);
  for (int trial = 0; trial < trials; ++trial) {
    const Constraints<P, M, N> constraints =
        RandomConstraintsGeneric<P, M, N>(rng);
    const int query = search_dist(rng);
    const int store = store_dist(rng);
    Constraints<P, M, N> transformed;
    for (const auto &r : constraints) {
      transformed.push_back(g.query.Apply(query, g.store.Apply(store, r)));
    }
    const int td = GaussJordanRREF<P, M, N>(&transformed);
    EXPECT_EQ(td, static_cast<int>(constraints.size()));
    const T t0 = ApplyConstraintsToTensor<P, M, N, N, N>(constraints, tensor);
    const T t1 = ApplyConstraintsToTensor<P, M, N, N, N>(transformed, tensor);
    EXPECT_TRUE(AreTensorsIsomorphic(t0, t1))
        << "P=" << P << " M=" << M << " N=" << N << " trial=" << trial
        << " query=" << query << " store=" << store;
  }
}

// q = 4 (P=2, M=2) and q = 9 (P=3, M=2): the base field is now a proper
// extension, so the query side is the full Gal(𝔽_{q^N}/𝔽_P) of size M·N (with
// genuinely semilinear elements) and the store side is the projective scaling
// quotient of size (q^N−1)/(q−1).
TEST(SymmetryGroupTest, GroupActionIsValidOverFq) {
  CheckGroupGeneric<2, 2, 2>();
  CheckGroupGeneric<2, 2, 3>();
  CheckGroupGeneric<3, 2, 2>();
}

TEST(SymmetryGroupTest, QueryStoreFactorizationIsUniqueOverFq) {
  CheckFactorizationGeneric<2, 2, 2>();
  CheckFactorizationGeneric<2, 2, 3>();
  CheckFactorizationGeneric<3, 2, 2>();
}

TEST(SymmetryGroupTest, ConstraintCommutesWithSymmetryOverFq) {
  std::mt19937_64 rng;
  CheckConstraintSymmetryIsoGeneric<2, 2, 2>(rng, 30);
  CheckConstraintSymmetryIsoGeneric<3, 2, 2>(rng, 12);
}

} // namespace
} // namespace extfield
