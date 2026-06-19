#pragma once

// A-side symmetry group for full (non-cyclic) polynomial multiplication over
// 𝔽_q with q = P^M.
//
//   T_full(N) : 𝔽_q[x]_<N × 𝔽_q[x]_<N → 𝔽_q[x]_<2N−1
//   T_full = Σ_{i,j < N} a_i ⊗ b_j ⊗ c_{i+j}.
//
// Automorphisms of T_full act on the A-side coefficient space 𝔽_q^N as the
// representation of GL(2, 𝔽_q) on degree-(N−1) binary forms F(X, Y) with
// p(x) = F(x, 1):
//   - substitution  S_c : x ↦ x+c for c ∈ 𝔽_q. S_c has order p in
//     characteristic p, and {S_c : c ∈ 𝔽_q} is the additive group (𝔽_q, +)
//     of order q. To generate it from a small set we use the M shifts
//     S_{y^k}, k = 0, …, M−1, where y^k is the F_q index P^k (the natural
//     𝔽_p-basis of 𝔽_q). For M = 1 this collapses to the single shift
//     S_1 = (x ↦ x+1).
//   - reverse  σ : a_i ↦ a_{N−1−i} (the x ↦ 1/x reversal), an involution.
//   - scaling  D_g : x ↦ g·x for g a generator of 𝔽_q^* (cyclic of order
//     q−1). Trivial for q = 2; nontrivial otherwise.
//
// Plus an extra *semilinear* symmetry whenever M ≥ 2: the Frobenius
// φ: a ↦ a^p generates Gal(𝔽_q/𝔽_p), cyclic of order M. T_full's structure
// constants live in {0, 1} ⊆ 𝔽_p, so the diagonal action (φ, φ, φ) on
// (A, B, C) fixes T_full. Combined with PGL we get
//
//     PΓL(2, 𝔽_q) = PGL(2, 𝔽_q) ⋊ Gal(𝔽_q/𝔽_p),
//     |PΓL(2, 𝔽_q)| = M · q(q²−1).
//
// We only need the action on constraint *subspaces* of the dual A-space, not
// on individual vectors, so the scalar centre {λ·I : λ ∈ 𝔽_q^*} of
// GL(2, 𝔽_q) acts trivially (it rescales every row uniformly, leaving the
// row span — and therefore every RREF key — unchanged). PGL is the projective
// quotient GL(2, 𝔽_q) / 𝔽_q^*, of order q(q²−1).
//
// Order of the working group PΓL(2, 𝔽_q) (for N ≥ 2; N = 1 collapses to the
// trivial group):
//   - q=2:    6.  q=3:   24.  q=5:    120.  q=7:    336.
//   - q=11: 1320. q=13: 2184.
//   - q=4:  120 (= 2·60).   q=8:  1512 (= 3·504).
//   - q=9: 1440 (= 2·720).  q=16: 16320 (= 4·4080).
//
// Factorisation: PGL goes on the *store* side; Gal goes on the *query* side.
// This is an *exact subgroup* factorisation — PGL is normal in PΓL with Gal a
// complement, so every (M, r) ∈ PΓL factors uniquely as (I, r) · (φ^{−r}(M), 0)
// = s⁻¹ · t with s⁻¹ ∈ Gal, t ∈ PGL — which is MORE than core/symmetry.h needs:
// the orbit map only requires the product set Query⁻¹ · Store to *cover* PΓL
// (a meet-in-the-middle), not a unique subgroup factorisation. Uniqueness is
// merely convenient here (one witness per orbit, no double-counting). The
// framework's "query after store" composition
// (`OrbitEnumeratorSlow::Transform`, `CanonicalFromWitness`) then reproduces
// the PΓL action (M, r)(v) = M · φ^r(v). Keeping PGL alone on the store side
// avoids an M× blow-up of `OrbitMap` memory; the cost is an M× slower
// `OrbitMap::Get` (M ≤ 4 in any q the repo cares about). For M = 1 Gal is
// trivial, query.Size() == 1, and the existing (P=2, M=1) BitVec hot path is
// untouched.
//
// Action on the dual A-space (store.Apply). A constraint is a functional in
// (𝔽_q^N)*. The transported constraint subspace under an automorphism with
// A-mode matrix M is the contragredient (M⁻¹)ᵀ. We build the group directly in
// this representation from the generators:
//   - G_{S_c} row i, col j = C(i, j) · c^(i−j) · (−1)^(i−j) for i ≥ j (else
//     0). For (P=2, M=1) and c=1 this collapses to [j ⊆ i] (XOR over
//     submasks).
//   - G_σ row i = e_{N−1−i} (σ is symmetric, so its dual is itself).
//   - G_D row i, col i = g⁻ⁱ (diagonal), where g generates 𝔽_q^*.
//
// Rows are stored as `GFVec<P, M, N>`, which transparently bit-packs the
// (P=2, M=1) row into a BitVec<N> (the legacy XOR/popcount hot path) and
// otherwise packs N base-field digits as `std::array<GF<P,M>, N>`. Mat is
// `std::array<GFVec<P,M,N>, N>`. BFS-closes from {S_{y^k} : k<M} ∪ {σ} ∪
// {D_g} (the scaling is omitted for q = 2 where 𝔽_q^* is trivial). Every
// product is canonicalised to the unique PGL representative whose first
// non-zero entry (row-major) is 1, so matrix-equality dedupe in IndexOf is
// exactly PGL-coset dedupe. For (P=2, M=1) the canonicalisation is a no-op
// (𝔽_2^* = {1}).
//
// Consequence of working in PGL: store.ApplyInverse(t, store.Apply(t, v))
// equals v only as a 1-D subspace, not as a Vec — Apply on the canonical
// representative produces some scalar multiple of the literal pre-image.
// Downstream code (OrbitMap.Set/Get, the verifier's RankMap) re-runs
// GaussJordanRREF after every Apply, which collapses the scalar back to a
// canonical key, so this projective composition is exactly what's wanted.

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/bit_vec.h"
#include "core/constraints.h"
#include "core/gf.h"
#include "core/gf_vec.h"

namespace full {

template <int P, int M, int N> class SymmetryGroup {
  static_assert(N >= 1);
  static_assert(P >= 2);
  static_assert(M >= 1);

public:
  static constexpr int kM = M;
  using Vec = GFVec<P, M, N>;

  // Query side: Gal(𝔽_q/𝔽_p), cyclic of order M, generated by the Frobenius
  // φ: a ↦ a^p. Element r ∈ {0, …, M−1} represents φ^r, acting entrywise on a
  // GFVec as v_i ↦ v_i^(p^r). Trivial for M = 1 (Gal(𝔽_p/𝔽_p) = {e}).
  struct QuerySet {
    using Elem = int;

    int Size() const { return M; }
    Elem At(int i) const { return i; }
    Elem Identity() const { return 0; }
    Vec ApplyInverse(Elem r, Vec v) const { return Apply((M - r) % M, v); }
    Vec Apply(Elem r, Vec v) const {
      if constexpr (M == 1) {
        return v;
      } else {
        const int exponent = static_cast<int>(IntPow(P, r));
        Vec out{};
        for (int i = 0; i < N; ++i) {
          out.Set(i, GF<P, M>::Pow(v[i], exponent));
        }
        return out;
      }
    }
  };

  // Store side: PGL(2, 𝔽_q) acting via Sym^(N−1), enumerated by BFS closure.
  class StoreSet {
  public:
    using Elem = int;

    StoreSet() {
      // Generators in the contragredient (dual) representation.
      std::vector<Mat> gens;
      // M translations x ↦ x + y^k jointly generate (𝔽_q, +).
      for (int k = 0; k < M; ++k) {
        gens.push_back(SubstitutionDual(BasisElement(k)));
      }
      gens.push_back(ReverseDual());
      // 𝔽_q^* is trivial for q = 2, so omit the scaling generator there.
      if constexpr (!(P == 2 && M == 1)) {
        gens.push_back(ScalingDual());
      }

      // Closure under right-multiplication by each generator, identity first.
      // Each product is canonicalised to the unique PGL representative (first
      // non-zero entry = 1), so matrix-equality dedupe in IndexOf is exactly
      // PGL-coset dedupe. For (P=2, M=1) the canonicalisation is a no-op.
      elems_.push_back(IdentityMat());
      for (std::size_t front = 0; front < elems_.size(); ++front) {
        for (const Mat &gen : gens) {
          const Mat product = Canonicalize(MatMul(elems_[front], gen));
          if (IndexOf(product) < 0) {
            elems_.push_back(product);
          }
        }
      }
      if constexpr (N == 1) {
        // For N = 1 every group element collapses to the identity (PGL acts
        // trivially on a 1-D coefficient space).
        CHECK_EQ(elems_.size(), 1);
      } else {
        constexpr std::size_t kQ = Ops::kQ;
        CHECK_EQ(elems_.size(), kQ * (kQ * kQ - 1))
            << "P=" << P << ", M=" << M << ", N=" << N;
      }

      // Precompute inverses: inverse_index_[i] is the j whose product
      // elems_[j]·elems_[i] is a scalar multiple of I (i.e., the identity in
      // PGL). Canonicalize collapses any λI to I, so we compare against
      // IdentityMat().
      const Mat id = IdentityMat();
      inverse_index_.assign(elems_.size(), -1);
      for (std::size_t i = 0; i < elems_.size(); ++i) {
        for (std::size_t j = 0; j < elems_.size(); ++j) {
          if (Canonicalize(MatMul(elems_[j], elems_[i])) == id) {
            inverse_index_[i] = static_cast<int>(j);
            break;
          }
        }
      }
    }

    int Size() const { return static_cast<int>(elems_.size()); }
    Elem At(int j) const { return j; }
    Elem Identity() const { return 0; }
    Vec Apply(Elem t, Vec v) const { return MatVec(elems_[t], v); }
    // Inverse action via the precomputed inverse element (this store is a
    // closed group, so the inverse is itself an element — see core/symmetry.h
    // on why that is convenient but not required).
    Vec ApplyInverse(Elem t, Vec v) const {
      return MatVec(elems_[inverse_index_[t]], v);
    }

  private:
    // Row representation: a single GFVec<P, M, N>. For (P=2, M=1) the
    // GFVec specialization bit-packs into a BitVec<N> (the XOR/popcount hot
    // path); otherwise it packs N base-field digits. Mat is N rows of Row.
    using Row = Vec;
    using Mat = std::array<Row, N>;
    using Ops = GF<P, M>;

    static Mat IdentityMat() {
      Mat m{};
      for (int i = 0; i < N; ++i) {
        m[i].Set(i, Ops::One());
      }
      return m;
    }

    // Pascal's triangle mod P, computed once per call (cheap for N ≤ 32).
    static std::array<std::array<uint8_t, N>, N> BinomTable() {
      std::array<std::array<uint8_t, N>, N> c{};
      for (int i = 0; i < N; ++i) {
        c[i][0] = 1;
        for (int j = 1; j <= i; ++j) {
          const int sum = static_cast<int>(c[i - 1][j - 1]) +
                          (j <= i - 1 ? static_cast<int>(c[i - 1][j]) : 0);
          c[i][j] = static_cast<uint8_t>(sum % P);
        }
      }
      return c;
    }

    // The F_q element y^k (the F_q index P^k), used as the k-th 𝔽_p-basis
    // shift generator for the additive group (𝔽_q, +).
    static Ops BasisElement(int k) {
      return Ops{static_cast<uint8_t>(IntPow(P, k))};
    }

    // Contragredient of x ↦ x + c: row i, col j = C(i, j) · c^(i−j) ·
    // (−1)^(i−j) for i ≥ j. For (P=2, M=1) and c = 1 this collapses to [j ⊆ i].
    static Mat SubstitutionDual(Ops c) {
      Mat m{};
      if constexpr (P == 2 && M == 1) {
        (void)c; // 𝔽_2 has only c = 1; the [j ⊆ i] body is c-independent.
        using BV = BitVec<N>;
        for (int i = 0; i < N; ++i) {
          BV row = 0;
          for (int j = 0; j <= i; ++j) {
            if ((i & j) == j) {
              row = static_cast<BV>(row | static_cast<BV>(BV{1} << j));
            }
          }
          m[i] = Row{row};
        }
      } else {
        const auto bin = BinomTable();
        for (int i = 0; i < N; ++i) {
          Ops cpow = Ops::One(); // c^0 for the k = 0 (j = i) cell
          for (int k = 0; k <= i; ++k) {
            const int j = i - k;
            Ops v{static_cast<uint8_t>(bin[i][j])};
            v = Ops::Mul(v, cpow);
            if (k & 1) {
              v = Ops::Neg(v);
            }
            m[i].Set(j, v);
            cpow = Ops::Mul(cpow, c);
          }
        }
      }
      return m;
    }

    // Contragredient of σ (involution σ = σᵀ): row i = e_{N−1−i}.
    static Mat ReverseDual() {
      Mat m{};
      for (int i = 0; i < N; ++i) {
        m[i].Set(N - 1 - i, Ops::One());
      }
      return m;
    }

    // A generator of 𝔽_q^* (cyclic of order q − 1). For q ≤ 16 the brute-force
    // scan is trivial. Returns Ops::One() for q = 2 (where 𝔽_q^* is trivial).
    static Ops FindFqStarGenerator() {
      constexpr int kQ = Ops::kQ;
      if constexpr (kQ == 2) {
        return Ops::One();
      } else {
        for (int g_idx = 2; g_idx < kQ; ++g_idx) {
          const Ops g{static_cast<uint8_t>(g_idx)};
          Ops cur = Ops::One();
          bool generates = true;
          for (int e = 0; e < kQ - 2; ++e) {
            cur = Ops::Mul(cur, g);
            if (cur == Ops::One()) {
              generates = false;
              break;
            }
          }
          if (generates) {
            return g;
          }
        }
        return Ops::Zero(); // unreachable for prime-power q
      }
    }

    // Contragredient of D_g (x ↦ g·x): diagonal with entry [i][i] = g⁻ⁱ.
    // The constructor omits this generator for (P=2, M=1) where the resulting
    // matrix would be the identity (𝔽_2^* is trivial); the body below
    // remains valid in that case (it returns IdentityMat).
    static Mat ScalingDual() {
      Mat m{};
      const Ops g_inv = Ops::Inverse(FindFqStarGenerator());
      Ops cur = Ops::One();
      for (int i = 0; i < N; ++i) {
        m[i].Set(i, cur);
        cur = Ops::Mul(cur, g_inv);
      }
      return m;
    }

    // (M v)_i = Σ_j M[i][j] · v_j. For (P=2, M=1) we use bitwise AND +
    // popcount parity; otherwise per-cell GF arithmetic.
    static Vec MatVec(const Mat &m, Vec v) {
      if constexpr (P == 2 && M == 1) {
        using BV = BitVec<N>;
        BV out = 0;
        const BV in = v.data;
        for (int i = 0; i < N; ++i) {
          if (__builtin_parityll(
                  static_cast<unsigned long long>(m[i].data & in))) {
            out = static_cast<BV>(out | static_cast<BV>(BV{1} << i));
          }
        }
        return Vec{out};
      } else {
        Vec out{};
        for (int i = 0; i < N; ++i) {
          Ops sum{};
          for (int j = 0; j < N; ++j) {
            sum = Ops::Add(sum, Ops::Mul(m[i][j], v[j]));
          }
          out.Set(i, sum);
        }
        return out;
      }
    }

    // Canonicalise an invertible matrix to the unique scalar representative
    // whose first non-zero entry (row-major) is 1. This is the canonical form
    // of m's coset in PGL(2, 𝔽_q) = GL(2, 𝔽_q) / 𝔽_q^*. For (P=2, M=1)
    // this is a no-op (𝔽_2^* = {1}). Invertibility guarantees row 0 has a
    // non-zero entry, so the pivot search succeeds.
    static Mat Canonicalize(Mat m) {
      if constexpr (P == 2 && M == 1) {
        return m;
      } else {
        int piv = 0;
        while (m[0][piv].value == 0) {
          ++piv;
          CHECK_LT(piv, N);
        }
        const Ops pivot = m[0][piv];
        if (pivot == Ops::One())
          return m;
        const Ops inv = Ops::Inverse(pivot);
        for (int i = 0; i < N; ++i) {
          m[i] *= inv;
        }
        return m;
      }
    }

    // (A·B)[i][k] = Σ_j A[i][j] · B[j][k]. Equivalently row i of A·B equals
    // Σ_j A[i][j] · (row j of B), which on the non-F₂ path uses the GFVec
    // scalar-multiply / vector-add operators directly.
    static Mat MatMul(const Mat &a, const Mat &b) {
      Mat out{};
      if constexpr (P == 2 && M == 1) {
        using BV = BitVec<N>;
        for (int i = 0; i < N; ++i) {
          BV row = 0;
          for (int j = 0; j < N; ++j) {
            if ((a[i].data >> j) & 1) {
              row = static_cast<BV>(row ^ b[j].data);
            }
          }
          out[i] = Row{row};
        }
      } else {
        for (int i = 0; i < N; ++i) {
          Row row{};
          for (int j = 0; j < N; ++j) {
            row += a[i][j] * b[j];
          }
          out[i] = row;
        }
      }
      return out;
    }

    int IndexOf(const Mat &m) const {
      for (std::size_t i = 0; i < elems_.size(); ++i) {
        if (elems_[i] == m) {
          return static_cast<int>(i);
        }
      }
      return -1;
    }

    std::vector<Mat> elems_;         // distinct group matrices, elems_[0] == I
    std::vector<int> inverse_index_; // element index -> index of its inverse
  };

  QuerySet query;
  StoreSet store;
};

} // namespace full
