#pragma once

// A-side symmetry group for cyclic polynomial multiplication over 𝔽_q with
// q = P^M.
//
//   T_cyc(N) : 𝔽_q[x]/(xᴺ−1) × 𝔽_q[x]/(xᴺ−1) → 𝔽_q[x]/(xᴺ−1)
//   T_cyc = Σ_{i,j < N} a_i ⊗ b_j ⊗ c_{(i+j) mod N}.
//
// Write R := 𝔽_q[x]/(xᴺ−1). Three families of A-side tensor autos:
//   - multiplication  L_p : a ↦ p·a for p ∈ R^*. Paired with (L_{p⁻¹}, I)
//     on (B, C): (p·a)·(p⁻¹·b) = a·b. The familiar "shift by r" is the
//     special case p = xʳ; constraining to ⟨x⟩ ⊂ R^* throws away all units
//     not of the form xʳ (most of R^*: by CRT, R ≅ ∏_d 𝔽_{q^d} over the
//     orbits of ·q on ℤ/N, so R^* ≅ ∏_d 𝔽_{q^d}^*, generally much bigger
//     than ⟨x⟩ ≅ ℤ/N).
//   - ring automorphism  φ ∈ Aut_{𝔽_q}(R), applied diagonally (φ, φ, φ).
//     Each φ is determined by φ(x) = y ∈ R with yᴺ = 1 and the powers
//     1, y, y², …, y^(N-1) linearly independent over 𝔽_q (equivalently,
//     𝔽_q[y] = R). The subgroup x ↦ xᵘ for u ∈ (ℤ/N)^* is only the
//     diagonal piece; Aut(R) also permutes equal-dimension CRT factors and
//     applies an independent Frobenius inside each factor.
//   - Galois  σ ∈ Gal(𝔽_q/𝔽_p), applied diagonally. Semilinear, only
//     non-trivial when M ≥ 2. T_cyc's structure constants live in
//     {0, 1} ⊆ 𝔽_p so the diagonal (σ, σ, σ) action fixes T_cyc.
//
// Multiplications and ring autos form a semidirect product
//   R^* ⋊ Aut(R)
// inside the linear group of R. The scalar centre {λ·1 : λ ∈ 𝔽_q^*} ⊂ R^*
// acts as scalar matrices on R, so for constraint-subspace canonicalisation
// (rows ↔ row span) it's invisible. The full projective A-side group is
//
//   G := (R^* ⋊ Aut(R) ⋊ Gal) / 𝔽_q^*,
//   |G| = M · |R^*| · |Aut(R)| / |𝔽_q^*|.
//
// |Aut(R)| factors over the CRT decomposition R ≅ ∏ 𝔽_{q^{d_i}}^{e_i} (with
// multiplicities e_i for equal extension degrees) as
//   |Aut(R)| = ∏_i d_i^{e_i} · e_i!.
//
// The two sides of the meet-in-the-middle split (see below) are
//
//   Query := Aut(R) ⋊ Gal,   |Query| = M · |Aut(R)|,
//   Store  := R^* / 𝔽_q^*,     |Store|  = |R^*| / (q − 1).
//
// Concrete sizes for (P=2, M=1) (current AGL = N·φ(N) for comparison):
//
//   N | xᴺ−1 factor degs  | |Store|=|R^*| | |Query|=|Aut(R)| | AGL
//  ---+-------------------+---------------+-------------------+-----
//   1 |  -                |            1  |               1   |   1
//   2 |  1+1 (char-2)     |            2  |               1   |   2
//   3 |  1+2              |            3  |               2   |   6
//   5 |  1+4              |           15  |               4   |  20
//   7 |  1+3+3            |           49  |              18   |  42
//   9 |  1+2+6            |          189  |              12   |  54
//  11 |  1+10             |        1 023  |              10   | 110
//  13 |  1+12             |        4 095  |              12   | 156
//  15 |  1+2+4+4+4        |       10 125  |             768   | 120
//
// Larger q (sizes computed once; verified in symmetry_test.cc), as
// (|Store|, |Query|):
//   (P=2, M=2, N=1):  R = 𝔽_4         → (1, 2).
//   (P=2, M=2, N=2):  R = 𝔽_4[t]/t²   → (4, 6).
//   (P=2, M=2, N=3):  R ≅ 𝔽_4³        → (9, 12).
//   (P=2, M=3, N=2):  R = 𝔽_8[t]/t²   → (8, 21).
//   (P=3, M=1, N=2):  R ≅ 𝔽_3²        → (2, 2).
//   (P=3, M=2, N=2):  R ≅ 𝔽_9²        → (8, 4).
//   (P=5, M=1, N=2):  R ≅ 𝔽_5²        → (4, 2).
//
// Factorisation across Query⁻¹ · Store: the semilinear Galois piece and the
// ring autos Aut(R) go on Query; the multiplications R^*/𝔽_q^* go on Store.
// Query is the subgroup H = Aut(R) ⋊ Gal (Gal normalises Aut(R): conjugating
// an 𝔽_q-linear ring auto by the Frobenius gives another 𝔽_q-linear ring auto),
// Store is the normal subgroup N = R^*/𝔽_q^*, and H ∩ N = {e}. So G = H · N is
// an exact (Zappa–Szép) factorisation — more than the framework needs
// (core/symmetry.h only asks the product set Query⁻¹ · Store = H⁻¹ · N = H · N
// to cover G); the uniqueness just keeps one witness per orbit. As H is a
// subgroup its action set is inverse-closed, so the enumerator's Query⁻¹·Store
// and OrbitMap::Get's Query⁻¹·Store coincide and both cover G.
//
// Why this split (vs. putting Aut(R) on Store too): OrbitMap::Set materialises
// |Store| rows per orbit while OrbitMap::Get probes the map |Query| times.
// Moving Aut(R) onto Query shrinks OrbitMap memory by a factor of |Aut(R)| at
// the cost of an |Aut(R)|× slower Get — a genuine √|G| rebalancing rather than
// the old "one big factor on store, tiny factor on query" shape. At N = 15
// the store drops from ≈ 7.8M to |R^*| = 10 125. For (P=2, M=1) the Galois
// factor is trivial (the Query index reduces to the auto index and the
// Frobenius step is a no-op), so the existing BitVec hot path is undisturbed.
//
// Action on the dual A-space. A constraint is a functional in (𝔽_q^N)*. The
// transported constraint subspace under an automorphism with A-mode matrix M
// is the contragredient (M⁻¹)ᵀ. Both sides build their group in this dual
// representation:
//   1. brute-force enumerate R^* by iterating all qᴺ length-N coefficient
//      vectors p and testing rank-N of the circulant L_p (column j = p · xʲ).
//      Keep coset reps mod 𝔽_q^* — the highest-order non-zero coefficient is
//      1. For (P=2, M=1) the filter is a no-op (𝔽_2^* = {1}).
//   2. brute-force enumerate Aut(R) by iterating qᴺ candidate images y of x
//      and accepting those with yᴺ = 1 and {1, y, …, y^(N-1)} linearly
//      independent.
//   3. Store builds one dual matrix per canonical unit p — the contragredient
//      of L_p (identity auto, y = x) — and canonicalises it to the unique
//      scalar representative whose first non-zero entry (row-major) is 1.
//      Distinct units give distinct canonical duals, so no IndexOf scan is
//      needed; the inverse table is built by inverting each Mat with
//      Gauss-Jordan and looking it up in a boost::unordered_flat_map<Mat,int>.
//   4. Query builds, per ring auto y, the contragredient of φ_y (unit p = 1)
//      and its *literal* (un-canonicalised) matrix inverse, so that the
//      composite action below round-trips exactly. The Galois factor acts
//      entrywise.
//
// Query element encoding. A Query element index e ∈ [0, M·|Aut(R)|)
// decodes auto-major as a = e / M (the ring auto) and r = e % M (the Galois
// power φ^r). Its action on a dual functional v is
//   Apply(e, v)        = M_a · Frobenius^r(v),
//   ApplyInverse(e, v) = Frobenius^{(M−r) mod M}( M_a⁻¹ · v ),
// where M_a is the auto's literal dual matrix. Because M_a⁻¹ is the literal
// inverse, ApplyInverse exactly inverts Apply: Frob^{M−r}(M_a⁻¹·M_a·Frob^r(v))
// = Frob^M(v) = v. e = 0 is the identity (a = 0 is the identity auto y = x,
// r = 0). core/orbit_enumerator_slow.h relies on this exact inverse.
//
// Consequence of working in PGL-style on the Store side:
// store.ApplyInverse(t, store.Apply(t, v)) equals v only as a 1-D subspace,
// not as a Vec — Apply on the canonical representative produces some scalar
// multiple of the literal pre-image. Downstream code (OrbitMap.Set/Get, the
// verifier's RankMap) re-runs GaussJordanRREF after every Apply, which
// collapses the scalar back to a canonical key, so this projective
// composition is exactly what's wanted.
//
// Rows are stored as `GFVec<P, M, N>`, which bit-packs the (P=2, M=1)
// row into a BitVec<N> (the legacy XOR/popcount hot path) and otherwise
// packs N base-field digits. Mat is `std::array<GFVec<P,M,N>, N>`. MatVec
// uses BitVec AND + popcount for (P=2, M=1); otherwise per-cell GF arithmetic.

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

namespace cyclic {

template <int P, int M, int N> class SymmetryGroup {
  static_assert(N >= 1);
  static_assert(P >= 2);
  static_assert(M >= 1);

public:
  static constexpr int kM = M;
  using Vec = GFVec<P, M, N>;

private:
  // Shared linear algebra over R = 𝔽_q[x]/(xᴺ−1) and its dual A-space, used by
  // both Query (Aut(R) dual matrices) and Store (R^* dual matrices). Members
  // are public within Detail, but Detail itself is private to the group.
  struct Detail {
    // Row representation. For (P=2, M=1) the GFVec specialization bit-packs
    // into a BitVec<N>; otherwise N base-field digits. Mat is N rows of Row.
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

    // The monomial x ∈ R, as a coefficient Vec (the image of x under the
    // identity ring auto). For N = 1 this reduces to 1 (x ≡ 1 mod x−1).
    static Vec XVec() { return RightCyclicShift(OneVec(), 1); }

    // Cyclic convolution: out_k = Σ_{i+j ≡ k mod N} a_i · b_j. Equivalent to
    // polynomial multiplication mod (xᴺ−1).
    static Vec CyclicConv(Vec a, Vec b) {
      if constexpr (P == 2 && M == 1) {
        using BV = BitVec<N>;
        const BV a_data = a.data;
        const BV b_data = b.data;
        BV res = 0;
        for (int i = 0; i < N; ++i) {
          if (!((a_data >> i) & 1))
            continue;
          const BV shifted =
              (i == 0) ? b_data
                       : static_cast<BV>(((b_data << i) | (b_data >> (N - i))) &
                                         kBitVecAllOnes<N>);
          res = static_cast<BV>(res ^ shifted);
        }
        return Vec{res};
      } else {
        Vec out{};
        for (int i = 0; i < N; ++i) {
          if (a[i].value == 0)
            continue;
          for (int j = 0; j < N; ++j) {
            if (b[j].value == 0)
              continue;
            const int k = (i + j) % N;
            out.Set(k, Ops::Add(out[k], Ops::Mul(a[i], b[j])));
          }
        }
        return out;
      }
    }

    // Cyclic right-shift of p by j positions: out_i = p_{(i-j) mod N}.
    // Equivalent to multiplying p by xʲ in R = 𝔽_q[x]/(xᴺ−1).
    static Vec RightCyclicShift(Vec p, int j) {
      const int jj = ((j % N) + N) % N;
      if constexpr (P == 2 && M == 1) {
        using BV = BitVec<N>;
        const BV d = p.data;
        const BV res = (jj == 0)
                           ? d
                           : static_cast<BV>(((d << jj) | (d >> (N - jj))) &
                                             kBitVecAllOnes<N>);
        return Vec{res};
      } else {
        Vec out{};
        for (int i = 0; i < N; ++i) {
          out.Set(i, p[(i - jj + N) % N]);
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
          result = CyclicConv(result, base);
        }
        base = CyclicConv(base, base);
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
          const Ops c = a[row][col];
          if (c.value == 0)
            continue;
          // a[row] -= c · a[col]; inv[row] -= c · inv[col].
          a[row] -= c * a[col];
          inv[row] -= c * inv[col];
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
        for (int j = 0; j < N; ++j) {
          const Vec col = RightCyclicShift(p, j);
          for (int i = 0; i < N; ++i) {
            lp[i].Set(j, col[i]);
          }
        }
        if (IsInvertible(lp)) {
          units.push_back(p);
        }
      }
      return units;
    }

    // Enumerate the images φ(x) over all ring autos φ ∈ Aut_{𝔽_q}(R). y is
    // accepted iff yᴺ = 1 in R AND the matrix [y⁰ | y¹ | … | y^(N-1)] has
    // rank N (equivalently 𝔽_q[y] = R, so x ↦ y extends to a bijective
    // ring auto). The identity auto y = x is always accepted and, being the
    // lex-smallest non-constant image, is enumerated first (index 0).
    static std::vector<Vec> EnumerateAutoImages() {
      std::vector<Vec> autos;
      const Vec one_vec = OneVec();
      for (uint64_t code = 0; code < kQN; ++code) {
        const Vec y = DecodeGFVec<P, M, N>(code);
        if (y.IsZero())
          continue;
        if (!(PolyPow(y, N) == one_vec))
          continue;
        std::array<Vec, N> y_powers;
        y_powers[0] = one_vec;
        for (int j = 1; j < N; ++j) {
          y_powers[j] = CyclicConv(y_powers[j - 1], y);
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

    // Build the dual-A matrix for the element a ↦ p · φ_y(a) with φ_y: x ↦ y.
    // The primal matrix has column j equal to (p · yʲ); we construct its
    // transpose (rows are (p · yʲ)) directly, then invert. The result is
    //   M_dual = (M_primal⁻¹)ᵀ = (M_primalᵀ)⁻¹.
    static Mat BuildDualMatrix(Vec p, Vec y) {
      Mat primal_T{};
      Vec yj = OneVec();
      for (int j = 0; j < N; ++j) {
        primal_T[j] = CyclicConv(p, yj);
        if (j + 1 < N) {
          yj = CyclicConv(yj, y);
        }
      }
      return MatInverse(primal_T);
    }

    // Canonicalise an invertible matrix to the unique scalar representative
    // whose first non-zero entry (row-major) is 1. For (P=2, M=1) it's a
    // no-op. Invertibility guarantees row 0 has a non-zero entry, so the
    // pivot search succeeds.
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
        // Literal (un-canonicalised) dual matrix of φ_y (unit p = 1) and its
        // literal inverse, so the round trip below is exact.
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
      // Dual matrix of L_p (multiplication by the unit p, identity auto y = x).
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
      // Inverse table: each elems_[i]'s matrix inverse, canonicalised, looked
      // up via boost::unordered_flat_map. O(|elems_| · N³). R^*/𝔽_q^* is a
      // closed subgroup, so every inverse is itself an enumerated element.
      boost::unordered_flat_map<Mat, int, typename Detail::MatHash>
          matrix_to_index;
      for (std::size_t i = 0; i < elems_.size(); ++i) {
        matrix_to_index[elems_[i]] = static_cast<int>(i);
      }
      CHECK_EQ(matrix_to_index.size(), elems_.size())
          << "duplicate canonical matrices: P=" << P << " M=" << M
          << " N=" << N;

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
            << " i=" << i;
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

} // namespace cyclic
