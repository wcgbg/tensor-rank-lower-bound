#pragma once

// A-side symmetry group for extension-field multiplication 𝔽_{q^N} =
// 𝔽_q[x]/p(x) over the base field 𝔽_q with q = P^M (P prime, M ≥ 1).
//
//   T_ext(N) : 𝔽_{q^N} × 𝔽_{q^N} → 𝔽_{q^N},  (f, g) ↦ f · g mod p.
//
// Two families of A-side automorphisms act on the coefficient space 𝔽_q^N
// (identified with the field, basis 1, x, …, x^{N−1}):
//   - multiplicative scaling  ·α : a ↦ α · a, for α ∈ 𝔽_{q^N}^*. As a tensor
//     symmetry, (·α, ·β, ·αβ) preserves T_ext because (α a)(β b) = αβ (ab).
//     The scalings form a cyclic group of order q^N − 1.
//   - ring automorphism / Galois  σ : a ↦ a^P. A ring automorphism of 𝔽_{q^N}
//     (σ(ab) = σ(a)σ(b)), so the diagonal (σ, σ, σ) preserves T_ext, with
//     σ α σ⁻¹ = α^P. The full ring-automorphism group is Gal(𝔽_{q^N}/𝔽_P) =
//     ⟨a ↦ a^P⟩, cyclic of order M·N — NOT N. It splits as
//         Gal(𝔽_{q^N}/𝔽_P)   ⊃   Gal(𝔽_{q^N}/𝔽_q)   (the 𝔽_q-linear part)
//            order M·N                order N (= ⟨a ↦ a^q⟩)
//     with quotient Gal(𝔽_q/𝔽_P) (order M, the base-field Frobenius, acting
//     𝔽_q-semilinearly). This base-field Galois factor mirrors the third
//     "Galois" family of cyclic/symmetry.h and the extra semilinear symmetry
//     of full/symmetry.h; it is trivial at M = 1.
//
// Effective (projective) group. A constraint is a functional and an orbit key
// is a row span (column-reversed RREF), so the scalar centre 𝔽_q^* — the
// base-field scalars λ, acting as the scalar matrix λ·I (mult-by-α is scalar on
// 𝔽_q^N iff α ∈ 𝔽_q) — rescales every row uniformly and is invisible to
// canonicalisation. The kernel of the action on subspaces is exactly 𝔽_q^*, so
// quotienting it out gives
//
//   G_A = (𝔽_{q^N}^* / 𝔽_q^*) ⋊ Gal(𝔽_{q^N}/𝔽_P),
//   |G_A| = (q^N − 1)/(q − 1) · M · N,
//
// the same projective quotient cyclic (… / 𝔽_q^*) and full (PGL = GL / 𝔽_q^*)
// already take. At (P, M) = (2, 1) the quotient is trivial (𝔽_q^* = {1}) and
// |G_A| = (2^N − 1) · N.
//
// Completeness. The A-side autotopism equation h(ab) = f(a)g(b) forces
// f = (·c) ∘ H with c ∈ 𝔽_{q^N}^* and H a ring automorphism, so G_A above is
// exactly the autotopism projection — nothing larger preserves the whole
// tensor. The trace form makes 𝔽_{q^N} a Frobenius algebra: lowering the
// C-index by the Gram matrix G_{ij} = Tr(x^{i+j}) makes T totally symmetric,
// giving rank_A(U) = rank_C(G·U). That transpose duality adds NO new A-side
// symmetry — it only swaps A- and C-constraints and closes back to the
// identity (the naive candidate G⁻¹·M_c·G is not a symmetry; it spuriously
// generates all of GL(N, 𝔽_q), which would falsely merge genuinely distinct
// constraint classes). Empirically, for P = 2 (N ≤ 4) and P = 3 (N ≤ 3) the
// G_A orbits coincide one-for-one with the full constrained-tensor isomorphism
// classes (AreTensorsIsomorphic over GL³), so G_A is the complete equivalence.
//
// Query⁻¹ · Store factorization (this implementation). SEARCH = the
// field-automorphism group Gal(𝔽_{q^N}/𝔽_P) = ⟨a↦a^P⟩ (size M·N): element k is
// σ_k: a ↦ a^{P^k}. Its dual action on a functional w is D_k · Φ_k(w), where
// Φ_k entrywise-raises the 𝔽_q coordinates to the P^{k mod M} power and
// D_k = (A_k⁻¹)ᵀ with A_k's column i = (xⁱ)^{P^k}. At M = 1, Φ_k is the
// identity and σ_k reduces to the pure-matrix Frobenius dual (query size N,
// byte-identical to the legacy code). STORE = the projective scalings
// 𝔽_{q^N}^*/𝔽_q^* (size (q^N−1)/(q−1)): the coset reps whose leading
// 𝔽_q-coordinate is 1. These two subgroups give an exact factorisation —
// every group element writes *uniquely* as (scaling)·σ^k — and the framework's
// "query after store" composition (OrbitEnumeratorSlow::Transform) reproduces
// the full G_A action. Uniqueness is more than core/symmetry.h requires (it
// only needs the product set Query⁻¹ · Store to cover G_A, a
// meet-in-the-middle, not a subgroup factorisation); here it conveniently keeps
// one witness per orbit.
//
// NOT pure entrywise Frobenius. Because the lex-smallest irreducible p
// generally has coefficients in 𝔽_q ∖ 𝔽_P, the bare diagonal (Φ_k, Φ_k, Φ_k)
// does NOT fix T (it maps mult-mod-p to mult-mod-p^(P)); the matrix D_k is the
// 𝔽_q-linear correction that makes σ_k the genuine field automorphism. This is
// why extfield cannot copy the pure entrywise-Frobenius query side of
// cyclic/full (whose 𝔽_P-rational structure constants do admit it).
//
// Element encoding. QuerySet::Elem = int k ∈ [0, M·N) for σ_k (k = 0 ≡ id).
// StoreSet::Elem = int j indexing the j-th scaling coset rep in lex order of
// α (j = 0 ≡ α = 1, the identity). For (P, M) == (2, 1) every nonzero α is its
// own coset rep, so j ↔ α = j + 1 as a bit pattern (the legacy convention) and
// the query side is the legacy Frobenius dual — both unchanged.
//
// Reachable fields: GF<P, M> caps q = P^M ≤ 16, so M ≥ 2 means
// q ∈ {4, 8, 16} (P = 2) or q = 9 (P = 3); N ≤ 16, and the slow enumerator
// further caps q^N ≤ 2^30.

#include <array>
#include <cstdint>
#include <map>
#include <vector>

#include <ng-log/logging.h>

#include "core/bit_vec.h"
#include "core/gf.h"
#include "core/gf_vec.h"
#include "problems/extfield/irreducibles.h"

namespace extfield {

// --- F_2 fast path: legacy FieldMul/FieldInverse on BitVec<N> ---
// Kept as free functions: the (P, M) == (2, 1) field-arithmetic hot path reuses
// them, and problems/extfield/symmetry_test.cc exercises them directly.

template <int N>
BitVec<N> FieldMul(BitVec<N> a, BitVec<N> b, BitVec<N> modulus) {
  using BV = BitVec<N>;
  constexpr BV mask = kBitVecAllOnes<N>;
  BV result = 0;
  BV acc = static_cast<BV>(a & mask);
  for (int i = 0; i < N; ++i) {
    if ((b >> i) & 1) {
      result ^= acc;
    }
    const bool overflow = (acc >> (N - 1)) & 1;
    acc = static_cast<BV>((acc << 1) & mask);
    if (overflow) {
      acc ^= modulus;
    }
  }
  return static_cast<BV>(result & mask);
}

template <int N> BitVec<N> FieldInverse(BitVec<N> a, BitVec<N> modulus) {
  using BV = BitVec<N>;
  std::uint64_t e = (std::uint64_t{1} << N) - 2;
  BV result = 1;
  BV base = a;
  while (e) {
    if (e & 1) {
      result = FieldMul<N>(result, base, modulus);
    }
    base = FieldMul<N>(base, base, modulus);
    e >>= 1;
  }
  return result;
}

// --- 𝔽_q field + matrix arithmetic on GFVec<P, M, N> ---
//
// 𝔽_{q^N} elements (in the basis 1, x, …, x^{N−1}), matrix rows, and dual-A
// functionals all share the row type GFVec<P, M, N>. Matrices are N rows.
// For (P, M) == (2, 1) the GFVec specialization bit-packs a BitVec<N>, so the
// `if constexpr (P == 2 && M == 1)` branches collapse to the legacy
// XOR/AND/popcount/shift hot paths.
namespace extfield_internal {

template <int P, int M, int N> using Vec = GFVec<P, M, N>;
template <int P, int M, int N> using Mat = std::array<GFVec<P, M, N>, N>;

// x^i (i < N) as a field element.
template <int P, int M, int N> Vec<P, M, N> XPow(int i) {
  Vec<P, M, N> v{};
  v.Set(i, GF<P, M>::One());
  return v;
}

template <int P, int M, int N> Vec<P, M, N> OneVec() {
  return XPow<P, M, N>(0);
}

// Low coefficients of the lex-smallest monic irreducible of degree N over 𝔽_q
// (x^N implicit One()).
template <int P, int M, int N> Vec<P, M, N> ModulusVec() {
  Vec<P, M, N> v{};
  const auto coeffs = IrreduciblePolyCoeffs<P, M, N>();
  for (int i = 0; i < N; ++i) {
    v.Set(i, coeffs[i]);
  }
  return v;
}

// Multiply a field element by x mod p (one shift + reduction).
template <int P, int M, int N>
Vec<P, M, N> MulByX(Vec<P, M, N> a, const Vec<P, M, N> &p_low) {
  using Ops = GF<P, M>;
  if constexpr (P == 2 && M == 1) {
    using BV = BitVec<N>;
    constexpr BV mask = kBitVecAllOnes<N>;
    const BV d = a.data;
    const bool overflow = (d >> (N - 1)) & 1;
    BV r = static_cast<BV>((d << 1) & mask);
    if (overflow) {
      r = static_cast<BV>(r ^ p_low.data);
    }
    return Vec<2, 1, N>{r};
  } else {
    const Ops top = a[N - 1];
    Vec<P, M, N> r{};
    for (int k = N - 1; k >= 1; --k) {
      r.Set(k, a[k - 1]);
    }
    if (top.value != 0) {
      for (int k = 0; k < N; ++k) {
        r.Set(k, Ops::Sub(r[k], Ops::Mul(top, p_low[k])));
      }
    }
    return r;
  }
}

// Field multiplication mod p.
template <int P, int M, int N>
Vec<P, M, N> FieldMul(Vec<P, M, N> a, Vec<P, M, N> b,
                      const Vec<P, M, N> &p_low) {
  using Ops = GF<P, M>;
  if constexpr (P == 2 && M == 1) {
    return Vec<2, 1, N>{::extfield::FieldMul<N>(a.data, b.data, p_low.data)};
  } else {
    Vec<P, M, N> result{};
    Vec<P, M, N> acc = a; // a · x^0, already reduced
    for (int i = 0; i < N; ++i) {
      const Ops bi = b[i];
      if (bi.value != 0) {
        for (int k = 0; k < N; ++k) {
          result.Set(k, Ops::Add(result[k], Ops::Mul(bi, acc[k])));
        }
      }
      if (i + 1 < N) {
        acc = MulByX<P, M, N>(acc, p_low);
      }
    }
    return result;
  }
}

// a^exp in 𝔽_{q^N} via repeated squaring.
template <int P, int M, int N>
Vec<P, M, N> FieldPow(Vec<P, M, N> a, uint64_t exp, const Vec<P, M, N> &p_low) {
  Vec<P, M, N> result = OneVec<P, M, N>();
  Vec<P, M, N> base = a;
  while (exp > 0) {
    if (exp & 1) {
      result = FieldMul<P, M, N>(result, base, p_low);
    }
    base = FieldMul<P, M, N>(base, base, p_low);
    exp >>= 1;
  }
  return result;
}

// Multiplicative inverse in 𝔽_{q^N} (Fermat: a^{q^N − 2}). q^N ≤ 2^30 for the
// reachable configs (the slow-enumerator cap), so the exponent fits in u64.
template <int P, int M, int N>
Vec<P, M, N> FieldInverse(Vec<P, M, N> a, const Vec<P, M, N> &p_low) {
  if constexpr (P == 2 && M == 1) {
    return Vec<2, 1, N>{::extfield::FieldInverse<N>(a.data, p_low.data)};
  } else {
    constexpr uint64_t kQN = IntPow(GF<P, M>::kQ, N);
    const uint64_t exp = kQN - 2;
    return FieldPow<P, M, N>(a, exp, p_low);
  }
}

template <int P, int M, int N> Mat<P, M, N> IdentityMat() {
  Mat<P, M, N> m{};
  for (int i = 0; i < N; ++i) {
    m[i].Set(i, GF<P, M>::One());
  }
  return m;
}

// (M v)_i = Σ_j m[i][j] · v_j. (P, M) == (2, 1): AND + popcount parity per row.
template <int P, int M, int N>
Vec<P, M, N> MatVec(const Mat<P, M, N> &m, Vec<P, M, N> v) {
  using Ops = GF<P, M>;
  if constexpr (P == 2 && M == 1) {
    using BV = BitVec<N>;
    BV out = 0;
    const BV in = v.data;
    for (int i = 0; i < N; ++i) {
      if (__builtin_parityll(static_cast<unsigned long long>(m[i].data & in))) {
        out = static_cast<BV>(out | static_cast<BV>(BV{1} << i));
      }
    }
    return Vec<2, 1, N>{out};
  } else {
    Vec<P, M, N> out{};
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

// Transpose.
template <int P, int M, int N> Mat<P, M, N> Transpose(const Mat<P, M, N> &m) {
  Mat<P, M, N> t{};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      t[i].Set(j, m[j][i]);
    }
  }
  return t;
}

// Gauss-Jordan inverse over 𝔽_q on [A | I]. Writes A⁻¹ to *out and returns true
// on success; returns false if A is singular. (Mirrors cyclic's TryMatInverse.)
template <int P, int M, int N>
bool TryMatInverse(Mat<P, M, N> a, Mat<P, M, N> *out) {
  using Ops = GF<P, M>;
  Mat<P, M, N> inv = IdentityMat<P, M, N>();
  for (int col = 0; col < N; ++col) {
    int pivot = -1;
    for (int row = col; row < N; ++row) {
      if (a[row][col].value != 0) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) {
      return false;
    }
    if (pivot != col) {
      std::swap(a[col], a[pivot]);
      std::swap(inv[col], inv[pivot]);
    }
    const Ops piv = a[col][col];
    if (!(piv == Ops::One())) {
      const Ops inv_piv = Ops::Inverse(piv);
      a[col] *= inv_piv;
      inv[col] *= inv_piv;
    }
    for (int row = 0; row < N; ++row) {
      if (row == col) {
        continue;
      }
      const Ops c = a[row][col];
      if (c.value == 0) {
        continue;
      }
      a[row] -= c * a[col];
      inv[row] -= c * inv[col];
    }
  }
  *out = inv;
  return true;
}

// Dual of scaling by α: row i = α⁻¹ · x^i (the contragredient (M_α⁻¹)ᵀ).
template <int P, int M, int N>
Mat<P, M, N> ScalingDual(Vec<P, M, N> alpha, const Vec<P, M, N> &p_low) {
  Mat<P, M, N> m{};
  Vec<P, M, N> cur = FieldInverse<P, M, N>(alpha, p_low);
  for (int i = 0; i < N; ++i) {
    m[i] = cur;
    cur = MulByX<P, M, N>(cur, p_low);
  }
  return m;
}

// Scale a field element so its leading (highest-index) nonzero coordinate is 1
// — the canonical coset representative of α · 𝔽_q^*. For (P, M) == (2, 1) it's
// a no-op (the leading coordinate is already 1).
template <int P, int M, int N> Vec<P, M, N> NormalizeScaling(Vec<P, M, N> a) {
  const GF<P, M> lead = a.LeadingNonzero();
  if (lead.value == 0 || lead == GF<P, M>::One()) {
    return a;
  }
  return GF<P, M>::Inverse(lead) * a;
}

} // namespace extfield_internal

// --- The SymmetryGroup ---

template <int P, int M, int N> class SymmetryGroup {
  static_assert(N >= 1);
  static_assert(P >= 2 && P < 256);
  static_assert(M >= 1);

  using Vec_ = GFVec<P, M, N>;
  using Mat_ = extfield_internal::Mat<P, M, N>;

public:
  static constexpr int kM = M;
  using Vec = GFVec<P, M, N>;

  // Query side: the field-automorphism group Gal(𝔽_{q^N}/𝔽_P) = ⟨a↦a^P⟩, of
  // order M·N. Element k ∈ [0, M·N) is σ_k: a ↦ a^{P^k}; its dual action on a
  // functional w is D_k · Φ_k(w), where Φ_k entrywise-raises 𝔽_q coordinates to
  // the P^{k mod M} power and D_k = (A_k⁻¹)ᵀ with A_k's column i = (xⁱ)^{P^k}.
  // At M = 1, Φ_k is the identity and this is the pure-matrix Frobenius dual.
  class QuerySet {
  public:
    using Elem = int;

    QuerySet() {
      const Vec_ p_low = extfield_internal::ModulusVec<P, M, N>();
      for (int k = 0; k < M * N; ++k) {
        // A_k: column i = (x^i)^{P^k} (the 𝔽_q-linear part of σ_k).
        Mat_ a_k{};
        for (int i = 0; i < N; ++i) {
          const Vec_ col = extfield_internal::FieldPow<P, M, N>(
              extfield_internal::XPow<P, M, N>(i), IntPow(P, k), p_low);
          for (int r = 0; r < N; ++r) {
            a_k[r].Set(i, col[r]);
          }
        }
        Mat_ a_inv{};
        const bool ok = extfield_internal::TryMatInverse<P, M, N>(a_k, &a_inv);
        CHECK(ok) << "ExtField QuerySet: A_k singular (P=" << P << " M=" << M
                  << " N=" << N << " k=" << k << ")";
        dual_[k] = extfield_internal::Transpose<P, M, N>(a_inv);
      }
    }

    int Size() const { return M * N; }
    Elem At(int i) const { return i; }
    Elem Identity() const { return 0; }
    Vec ApplyInverse(Elem k, Vec v) const {
      return Apply((M * N - k) % (M * N), v);
    }
    Vec Apply(Elem k, Vec v) const {
      if constexpr (M == 1) {
        return extfield_internal::MatVec<P, M, N>(dual_[k], v);
      } else {
        // Entrywise Frobenius^{k mod M} on the 𝔽_q coordinates, then D_k.
        const int exp = static_cast<int>(IntPow(P, k % M));
        Vec tw{};
        for (int i = 0; i < N; ++i) {
          tw.Set(i, GF<P, M>::Pow(v[i], exp));
        }
        return extfield_internal::MatVec<P, M, N>(dual_[k], tw);
      }
    }

  private:
    std::array<Mat_, static_cast<std::size_t>(M *N)> dual_;
  };

  // Store side: projective scalings 𝔽_{q^N}^*/𝔽_q^*, the coset reps whose
  // leading 𝔽_q-coordinate is 1, of size (q^N − 1)/(q − 1).
  class StoreSet {
  public:
    using Elem = int;

    StoreSet() {
      const Vec_ p_low = extfield_internal::ModulusVec<P, M, N>();
      p_low_ = p_low;
      constexpr uint64_t kQ = GF<P, M>::kQ;
      constexpr uint64_t kQN = IntPow(kQ, N);
      for (uint64_t i = 1; i < kQN; ++i) {
        const Vec_ alpha = DecodeGFVec<P, M, N>(i);
        if (!(alpha.LeadingNonzero() == GF<P, M>::One())) {
          continue; // keep one representative per α · 𝔽_q^* coset
        }
        index_of_[EncodeGFVec<P, M, N>(alpha)] = static_cast<int>(reps_.size());
        reps_.push_back(alpha);
        duals_.push_back(extfield_internal::ScalingDual<P, M, N>(alpha, p_low));
      }
      CHECK(!reps_.empty());
      CHECK(reps_.size() * (kQ - 1) == kQN - 1);
    }

    int Size() const { return static_cast<int>(reps_.size()); }
    Elem At(int j) const { return j; }
    Elem Identity() const { return 0; }

    Vec Apply(Elem s, Vec v) const {
      return extfield_internal::MatVec<P, M, N>(duals_[s], v);
    }

    // Inverse action: the dual of scaling by reps_[s]⁻¹. This store is a closed
    // group, so that inverse scaling is itself an enumerated coset rep; we look
    // it up and apply its precomputed dual (see core/symmetry.h on why closure
    // is convenient but not required).
    Vec ApplyInverse(Elem s, Vec v) const {
      const Vec_ inv =
          extfield_internal::FieldInverse<P, M, N>(reps_[s], p_low_);
      const Vec_ rep = extfield_internal::NormalizeScaling<P, M, N>(inv);
      const auto it = index_of_.find(EncodeGFVec<P, M, N>(rep));
      CHECK(it != index_of_.end())
          << "ExtField StoreSet: inverse rep not found (P=" << P << " M=" << M
          << " N=" << N << " s=" << s << ")";
      return extfield_internal::MatVec<P, M, N>(duals_[it->second], v);
    }

  private:
    std::vector<Vec_> reps_;
    std::vector<Mat_> duals_;
    std::map<uint64_t, int> index_of_;
    Vec_ p_low_{};
  };

  QuerySet query;
  StoreSet store;
};

} // namespace extfield
