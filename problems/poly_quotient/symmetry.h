#pragma once

// A-side symmetry group for multiplication in the polynomial quotient ring
//   R := 𝔽_q[x]/(xᴺ − γ),   q = P^M,   γ ∈ 𝔽_p ⊆ 𝔽_q a scalar.
//
// This is the shared engine behind three problem families, which differ ONLY in
// the reduction constant γ (the value of xᴺ in R):
//   - cyclic      γ = 1   (xᴺ − 1)
//   - truncated   γ = 0   (xᴺ; x nilpotent — the short product)
//   - negacyclic  γ = −1  (xᴺ + 1; over characteristic 2, −1 = 1 = cyclic)
// It is parameterised by `kGamma`, the 𝔽_p index of γ (0, 1, or P−1). Because
// γ lies in the prime subfield, GFVec digit kGamma IS the element γ for every
// (P, M): index 0 → 0, 1 → 1, P−1 → −1.
//
//   T(N) : R × R → R,   (f, g) ↦ f · g mod (xᴺ − γ)
//   T = Σ_{i,j < N} a_i ⊗ b_j ⊗ c_k, where xⁱ⁺ʲ ≡ Σ_k (…) xᵏ mod (xᴺ − γ).
//
// Three families of A-side tensor autos act on coefficient space 𝔽_q^N, exactly
// as for cyclic (the construction never used xᴺ−1 specifically, only that R is a
// commutative 𝔽_q-algebra of the form 𝔽_q[x]/(monic)):
//   - multiplication  L_p : a ↦ p·a for p ∈ R^*. Paired with (L_{p⁻¹}, I) on
//     (B, C): (p·a)·(p⁻¹·b) = a·b. Valid for the truncated (non-reduced) ring
//     too: trunc(p·a·b) = trunc(p·trunc(a·b)) since the dropped tail has degree
//     ≥ N. The "shift by r" map is the special case p = xʳ.
//   - ring automorphism  φ ∈ Aut_{𝔽_q}(R), applied diagonally (φ, φ, φ). Each φ
//     is determined by φ(x) = y ∈ R with yᴺ = γ in R (the defining relation of
//     R) and 1, y, …, y^(N−1) linearly independent over 𝔽_q (𝔽_q[y] = R).
//   - Galois  σ ∈ Gal(𝔽_q/𝔽_p), applied diagonally. Semilinear, non-trivial
//     only for M ≥ 2. The structure constants live in 𝔽_p, so (σ, σ, σ) fixes T.
//
// Multiplications and ring autos form a semidirect product R^* ⋊ Aut(R) inside
// the linear group of R. The scalar centre {λ·1 : λ ∈ 𝔽_q^*} ⊂ R^* acts as
// scalar matrices, invisible to constraint-subspace canonicalisation. The full
// projective A-side group is
//
//   G := (R^* ⋊ Aut(R) ⋊ Gal) / 𝔽_q^*,
//   |G| = M · |R^*| · |Aut(R)| / |𝔽_q^*|.
//
// The two sides of the meet-in-the-middle split:
//   Query := Aut(R) ⋊ Gal,   |Query| = M · |Aut(R)|,
//   Store  := R^* / 𝔽_q^*,     |Store|  = |R^*| / (q − 1).
//
// Examples of the per-family unit/auto counts (P=2, M=1):
//   - cyclic (γ=1):     R^* = ∏ 𝔽_{q^d}^* over the q-cyclotomic cosets mod N.
//   - truncated (γ=0):  R = 𝔽_q[x]/(xᴺ), units are {p : p_0 ≠ 0}, so
//     |R^*| = (q−1)·q^(N−1) and |Store| = q^(N−1); autos are {y : y_0 = 0,
//     y_1 ≠ 0}, so |Aut(R)| = (q−1)·q^(N−2) (and 1 for N = 1, where x ≡ 0).
//   - negacyclic (γ=−1): R = 𝔽_q[x]/(xᴺ+1); the unit/auto counts follow the
//     factorisation of xᴺ+1 (cyclotomic cosets of the odd residues mod 2N for
//     p odd; identical to cyclic for p = 2 since xᴺ+1 = xᴺ−1).
//
// Factorisation across Query⁻¹ · Store: the Galois piece and Aut(R) go on Query
// (a subgroup, Gal normalises the 𝔽_q-linear Aut(R)); R^*/𝔽_q^* goes on Store
// (normal subgroup); H ∩ N = {e}. So G = H · N is an exact (Zappa–Szép)
// factorisation — more than core/symmetry.h needs (it asks only that the product
// set Query⁻¹ · Store cover G); the uniqueness just keeps one witness per orbit.
//
// Cost model: OrbitMap::Set materialises |Store| rows per orbit; OrbitMap::Get
// probes |Query| times. Moving Aut(R) onto Query rebalances toward √|G|. For
// (P=2, M=1) the Galois factor is trivial and the BitVec hot path is undisturbed.
//
// Action on the dual A-space. A constraint is a functional in (𝔽_q^N)*. Under an
// automorphism with A-mode matrix M the transported subspace is the
// contragredient (M⁻¹)ᵀ. Store builds, per canonical unit p, the contragredient
// of L_p and canonicalises it to the scalar rep whose first non-zero entry
// (row-major) is 1; Query builds, per ring auto y, the *literal* contragredient
// of φ_y and its literal inverse so the composite round-trips exactly. See the
// per-side comments below; this mirrors cyclic exactly — the only ring-specific
// pieces are Conv / MulByX / ReductionVec and the auto relation yᴺ = γ.
//
// Rows are `GFVec<P, M, N>` (bit-packed BitVec<N> for (P=2,M=1), else N base-
// field digits). Mat is `std::array<GFVec<P,M,N>, N>`.

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include <boost/unordered/unordered_flat_map.hpp>
#include <ng-log/logging.h>
#include <tbb/parallel_for.h>

#include "core/bit_vec.h"
#include "core/gf.h"
#include "core/gf_vec.h"

namespace poly_quotient {

// kGamma is the 𝔽_p index of the reduction constant γ (the value of xᴺ in R):
// 0 (truncated), 1 (cyclic), or P−1 (negacyclic, = −1). The three problem
// families alias this template with the appropriate kGamma.
template <int P, int M, int N, int kGamma> class SymmetryGroup {
  static_assert(N >= 1);
  static_assert(P >= 2);
  static_assert(M >= 1);
  static_assert(kGamma >= 0 && kGamma < P,
                "kGamma is an 𝔽_p index, so it must lie in [0, P)");
  // (kGamma is the 𝔽_p index of γ, the value of xᴺ in R: cyclic 1, truncated 0,
  //  negacyclic P−1.)

public:
  static constexpr int kM = M;
  using Vec = GFVec<P, M, N>;

private:
  // Shared linear algebra over R = 𝔽_q[x]/(xᴺ−γ) and its dual A-space, used by
  // both Query (Aut(R) dual matrices) and Store (R^* dual matrices). Members
  // are public within Detail, but Detail itself is private to the group.
  struct Detail {
    using Row = Vec;
    using Mat = std::array<Row, N>;
    using Ops = GF<P, M>;

    struct MatHash {
      size_t operator()(const Mat &m) const noexcept {
        const uint8_t *data = reinterpret_cast<const uint8_t *>(m.data());
        size_t bytes = sizeof(Mat);
        uint64_t h = bytes;
        while (bytes >= 8) {
          uint64_t chunk;
          std::memcpy(&chunk, data, 8);
          h = (chunk + (h << 6)) + (0x9e3779b9 + (h >> 2));
          data += 8;
          bytes -= 8;
        }
        if (bytes > 0) {
          uint64_t chunk = 0;
          std::memcpy(&chunk, data, bytes);
          h = (chunk + (h << 6)) + (0x9e3779b9 + (h >> 2));
        }
        return h;
      }
    };

    static constexpr uint64_t kQN = IntPow(Ops::kQ, N);

    // The reduction constant γ ∈ 𝔽_q (value of xᴺ in R), built from its 𝔽_p
    // index kGamma: index 0 → 0, 1 → 1, P−1 → −1, for every (P, M).
    static constexpr Ops Gamma() { return Ops{static_cast<uint8_t>(kGamma)}; }

    static Mat IdentityMat() {
      Mat m{};
      for (int i = 0; i < N; ++i) {
        m[i].Set(i, Ops::One());
      }
      return m;
    }

    // 1 ∈ R, as a coefficient Vec.
    static Vec OneVec() {
      Vec v{};
      v.Set(0, Ops::One());
      return v;
    }

    // The value of xᴺ in R, i.e. γ·1, as a coefficient Vec. An automorphism
    // x ↦ y must satisfy yᴺ = this in R.
    static Vec ReductionVec() {
      Vec v{};
      v.Set(0, Gamma());
      return v;
    }

    // Multiply p by x in R: out_0 = γ·p_{N−1}, out_i = p_{i−1} (i ≥ 1). For
    // N = 1 this is out_0 = γ·p_0 (x ≡ γ). The monomial x ∈ R is MulByX(1).
    static Vec MulByX(Vec p) {
      Vec out{};
      out.Set(0, Ops::Mul(Gamma(), p[N - 1]));
      for (int i = 1; i < N; ++i) {
        out.Set(i, p[i - 1]);
      }
      return out;
    }

    // The monomial x ∈ R, as a coefficient Vec.
    static Vec XVec() { return MulByX(OneVec()); }

    // Multiplication in R: schoolbook product to degree 2N−2, then fold
    // xᴺ⁺ᵗ ≡ γ·xᵗ, i.e. out_t = p_t + γ·p_{t+N} (and out_{N−1} = p_{N−1}).
    static Vec Conv(Vec a, Vec b) {
      if constexpr (P == 2 && M == 1) {
        using BV = BitVec<N>;
        const BV a_data = a.data;
        const BV b_data = b.data;
        BV res = 0;
        for (int i = 0; i < N; ++i) {
          if (!((a_data >> i) & 1))
            continue;
          const BV low = static_cast<BV>((b_data << i) & kBitVecAllOnes<N>);
          res = static_cast<BV>(res ^ low);
          if constexpr (kGamma == 1) {
            // Fold the part that overflowed past bit N−1 back to the low bits.
            // (kGamma == 0 drops it — the truncated product.)
            if (i != 0) {
              res = static_cast<BV>(res ^ static_cast<BV>(b_data >> (N - i)));
            }
          }
        }
        return Vec{res};
      } else {
        const Ops gamma = Gamma();
        std::array<Ops, 2 * N - 1> prod{};
        for (int i = 0; i < N; ++i) {
          if (a[i].value == 0)
            continue;
          for (int j = 0; j < N; ++j) {
            if (b[j].value == 0)
              continue;
            prod[i + j] = Ops::Add(prod[i + j], Ops::Mul(a[i], b[j]));
          }
        }
        Vec out{};
        for (int t = 0; t < N; ++t) {
          Ops val = prod[t];
          if (t + N <= 2 * N - 2) {
            val = Ops::Add(val, Ops::Mul(gamma, prod[t + N]));
          }
          out.Set(t, val);
        }
        return out;
      }
    }

    // y^exp in R via repeated squaring.
    static Vec PolyPow(Vec y, int exp) {
      Vec result = OneVec();
      Vec base = y;
      while (exp > 0) {
        if (exp & 1) {
          result = Conv(result, base);
        }
        base = Conv(base, base);
        exp >>= 1;
      }
      return result;
    }

    // Gauss-Jordan inverse over 𝔽_q on [A | I]. Returns true and writes A⁻¹
    // to *out on success; returns false if A is singular.
    static bool TryMatInverse(Mat a, Mat *out) {
      Mat inv = IdentityMat();
      for (int col = 0; col < N; ++col) {
        int pivot = -1;
        for (int row = col; row < N; ++row) {
          if (a[row][col].value != 0) {
            pivot = row;
            break;
          }
        }
        if (pivot < 0)
          return false;
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
          if (row == col)
            continue;
          const Ops co = a[row][col];
          if (co.value == 0)
            continue;
          a[row] -= co * a[col];
          inv[row] -= co * inv[col];
        }
      }
      *out = inv;
      return true;
    }

    static Mat MatInverse(const Mat &a) {
      Mat out;
      const bool ok = TryMatInverse(a, &out);
      CHECK(ok) << "MatInverse: singular matrix";
      return out;
    }

    static bool IsInvertible(const Mat &a) {
      Mat dummy;
      return TryMatInverse(a, &dummy);
    }

    // Enumerate units of R, kept as coset reps of R^*/𝔽_q^* (highest-order
    // non-zero coefficient = 1). For (P=2, M=1) the coset filter is a no-op.
    static std::vector<Vec> EnumerateCanonicalUnits() {
      std::vector<Vec> units;
      for (uint64_t code = 0; code < kQN; ++code) {
        const Vec p = DecodeGFVec<P, M, N>(code);
        if (p.IsZero())
          continue;
        if (!(p.LeadingNonzero() == Ops::One()))
          continue;
        Mat lp{};
        Vec col = p; // column j of L_p is p · xʲ.
        for (int j = 0; j < N; ++j) {
          for (int i = 0; i < N; ++i) {
            lp[i].Set(j, col[i]);
          }
          col = MulByX(col);
        }
        if (IsInvertible(lp)) {
          units.push_back(p);
        }
      }
      return units;
    }

    // Enumerate the images φ(x) over all ring autos φ ∈ Aut_{𝔽_q}(R). y is
    // accepted iff yᴺ = γ in R AND [y⁰ | y¹ | … | y^(N-1)] has rank N
    // (equivalently 𝔽_q[y] = R). The zero vector is NOT pre-filtered: for the
    // truncated ring at N = 1 the identity auto is x ↦ 0 (x ≡ 0 in 𝔽_q), which
    // satisfies y¹ = 0 = γ and {1} independent; for N ≥ 2, or γ ≠ 0, y = 0 fails
    // one of the two checks and is rejected. The identity auto y = x is always
    // accepted and enumerated as one of the images.
    static std::vector<Vec> EnumerateAutoImages() {
      std::vector<Vec> autos;
      const Vec reduction = ReductionVec();
      const Vec one_vec = OneVec();
      for (uint64_t code = 0; code < kQN; ++code) {
        const Vec y = DecodeGFVec<P, M, N>(code);
        if (!(PolyPow(y, N) == reduction))
          continue;
        std::array<Vec, N> y_powers;
        y_powers[0] = one_vec;
        for (int j = 1; j < N; ++j) {
          y_powers[j] = Conv(y_powers[j - 1], y);
        }
        Mat my{};
        for (int j = 0; j < N; ++j) {
          for (int i = 0; i < N; ++i) {
            my[i].Set(j, y_powers[j][i]);
          }
        }
        if (IsInvertible(my)) {
          autos.push_back(y);
        }
      }
      return autos;
    }

    // Build the dual-A matrix for a ↦ p · φ_y(a) with φ_y: x ↦ y. The primal
    // matrix has column j equal to (p · yʲ); we build its transpose (rows are
    // (p · yʲ)) and invert: M_dual = (M_primal⁻¹)ᵀ = (M_primalᵀ)⁻¹.
    static Mat BuildDualMatrix(Vec p, Vec y) {
      Mat primal_T{};
      Vec yj = OneVec();
      for (int j = 0; j < N; ++j) {
        primal_T[j] = Conv(p, yj);
        if (j + 1 < N) {
          yj = Conv(yj, y);
        }
      }
      return MatInverse(primal_T);
    }

    // Canonicalise an invertible matrix to the scalar rep whose first non-zero
    // entry (row-major) is 1. For (P=2, M=1) it's a no-op.
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

    // (M v)_i = Σ_j M[i][j] · v_j. (P=2, M=1) hot path: AND + popcount per row.
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

    // Entrywise Frobenius φ^r: v_i ↦ v_i^(P^r). No-op for M = 1.
    static Vec Frobenius(int r, Vec v) {
      if constexpr (M == 1) {
        return v;
      } else {
        const int exponent = static_cast<int>(IntPow(P, r));
        Vec out{};
        for (int i = 0; i < N; ++i) {
          out.Set(i, Ops::Pow(v[i], exponent));
        }
        return out;
      }
    }
  };

public:
  // Query side: Aut(R) ⋊ Gal. Element e ∈ [0, M·|Aut(R)|) decodes auto-major
  // as a = e / M (ring auto) and r = e % M (Galois power φ^r). Its action is
  //   Apply(e, v)        = M_a · Frobenius^r(v),
  //   ApplyInverse(e, v) = Frobenius^{(M−r) mod M}( M_a⁻¹ · v ),
  // with M_a the auto's literal (un-canonicalised) dual matrix and M_a⁻¹ its
  // literal inverse, so ApplyInverse exactly inverts Apply. For M = 1 the
  // Frobenius step is a no-op and e is just the auto index.
  class QuerySet {
  public:
    using Elem = int;
    using Mat = typename Detail::Mat;

    QuerySet() {
      const std::vector<Vec> auto_images = Detail::EnumerateAutoImages();
      num_autos_ = static_cast<int>(auto_images.size());
      const Vec one = Detail::OneVec();
      auto_dual_.resize(num_autos_);
      auto_dual_inv_.resize(num_autos_);
      for (int a = 0; a < num_autos_; ++a) {
        const Mat dual = Detail::BuildDualMatrix(one, auto_images[a]);
        auto_dual_[a] = dual;
        auto_dual_inv_[a] = Detail::MatInverse(dual);
      }
      if constexpr (N == 1) {
        CHECK_EQ(num_autos_, 1);
      }
    }

    int Size() const { return M * num_autos_; }
    Elem At(int i) const { return i; }
    Elem Identity() const { return 0; }

    Vec Apply(Elem e, Vec v) const {
      const int a = e / M;
      const int r = e % M;
      return Detail::MatVec(auto_dual_[a], Detail::Frobenius(r, v));
    }
    Vec ApplyInverse(Elem e, Vec v) const {
      const int a = e / M;
      const int r = e % M;
      return Detail::Frobenius((M - r) % M,
                               Detail::MatVec(auto_dual_inv_[a], v));
    }

  private:
    int num_autos_ = 0;
    std::vector<Mat> auto_dual_;     // dual matrix of each ring auto (literal)
    std::vector<Mat> auto_dual_inv_; // its literal inverse
  };

  // Store side: R^* / 𝔽_q^* (multiplications) acting on the dual A-space.
  class StoreSet {
  public:
    using Elem = int;
    using Mat = typename Detail::Mat;

    StoreSet() {
      const std::vector<Vec> canon_units = Detail::EnumerateCanonicalUnits();
      const std::size_t num_units = canon_units.size();
      elems_.resize(num_units);
      if (num_units > 1'000'000) {
        LOG(INFO) << "SymmetryGroup::StoreSet::StoreSet: num_elems="
                  << num_units;
        LOG(INFO) << "SymmetryGroup::StoreSet::StoreSet: Building "
                     "canonical matrices...";
      }
      const Vec x = Detail::XVec();
      tbb::parallel_for(std::size_t(0), num_units, [&](std::size_t i) {
        Mat dual = Detail::BuildDualMatrix(canon_units[i], x);
        dual = Detail::Canonicalize(dual);
        elems_[i] = dual;
      });

      if constexpr (N == 1) {
        CHECK_EQ(elems_.size(), 1u);
      }

      if (num_units > 1'000'000) {
        LOG(INFO) << "SymmetryGroup::StoreSet::StoreSet: Building "
                     "inverse table...";
      }
      boost::unordered_flat_map<Mat, int, typename Detail::MatHash>
          matrix_to_index;
      for (std::size_t i = 0; i < elems_.size(); ++i) {
        matrix_to_index[elems_[i]] = static_cast<int>(i);
      }
      CHECK_EQ(matrix_to_index.size(), elems_.size())
          << "duplicate canonical matrices: P=" << P << " M=" << M
          << " N=" << N << " kGamma=" << kGamma;

      if (num_units > 1'000'000) {
        LOG(INFO) << "SymmetryGroup::StoreSet::StoreSet: Computing "
                     "inverse indices...";
      }
      inverse_index_.assign(elems_.size(), -1);
      tbb::parallel_for(std::size_t(0), elems_.size(), [&](std::size_t i) {
        Mat inv = Detail::MatInverse(elems_[i]);
        inv = Detail::Canonicalize(inv);
        const auto it = matrix_to_index.find(inv);
        CHECK(it != matrix_to_index.end())
            << "inverse not found: P=" << P << " M=" << M << " N=" << N
            << " kGamma=" << kGamma << " i=" << i;
        inverse_index_[i] = it->second;
      });
      if (num_units > 1'000'000) {
        LOG(INFO) << "SymmetryGroup::StoreSet::StoreSet: Construction "
                     "complete.";
      }
    }

    int Size() const { return static_cast<int>(elems_.size()); }
    Elem At(int j) const { return j; }
    Elem Identity() const { return 0; }
    Vec Apply(Elem t, Vec v) const { return Detail::MatVec(elems_[t], v); }
    // Inverse action via the precomputed inverse element (this store is a
    // closed group, so the inverse is itself an element — see core/symmetry.h
    // on why that is convenient but not required).
    Vec ApplyInverse(Elem t, Vec v) const {
      return Detail::MatVec(elems_[inverse_index_[t]], v);
    }

  private:
    std::vector<Mat> elems_;         // canonical matrices, elems_[0] == I
    std::vector<int> inverse_index_; // element index -> index of its inverse
  };

  QuerySet query;
  StoreSet store;
};

} // namespace poly_quotient
