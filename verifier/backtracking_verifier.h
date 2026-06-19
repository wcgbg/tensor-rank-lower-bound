#pragma once

// Replay-verifier for the Substitution-with-Backtracking proof. Port of the
// matrix-mult proof_verifier/rank_lower_bound_backtracking_verifier.h.
//
// The prover (core/rank_lower_bound_backtracking.h) walks a deterministic DFS
// over the "minimal" extra constraints of the base subspace, emitting one
// BacktrackingProof record (in DFS pre-order) at each leaf where it reached the
// target bound. This verifier rebuilds the identical minimal-constraint list
// and re-walks the same DFS, consuming records in order: a node whose depth
// equals the next record's dfs_constraints_size is a recorded leaf (consume
// and return); otherwise it descends into all children. At each leaf it
// reconstructs the extended subspace from the record's mask, recovers the
// canonical orbit representative via the (query_elem, store_elem) witness
// (CanonicalFromWitness), looks its bound up in the plain canonical→rank map,
// and checks the substitution bound popcount(mask) + rank meets the claimed
// lower bound.

#include <bit>
#include <cstdint>
#include <vector>

#include <boost/unordered/unordered_flat_map.hpp>
#include <ng-log/logging.h>

#include "core/backtracking_proof.h"
#include "core/bit_vec.h"
#include "core/constraints.h"
#include "core/gf.h"
#include "core/gf_vec.h"

// Plain map from an orbit's canonical representative to its rank lower bound.
// The verifier keeps one entry per orbit, unlike the prover's square-root
// OrbitMap which materializes every Store-image.
template <class Problem>
using RankMap = boost::unordered_flat_map<
    Constraints<Problem::kP, Problem::kM, Problem::kNA>, int,
    ConstraintsHash<Problem::kP, Problem::kM, Problem::kNA>>;

// Recover the canonical representative c of query q's orbit from the recorded
// witness, using the OrbitMap hit equation
//   query.Apply(query_elem, q) ≡ store.Apply(store_elem, c)   (as subspaces)
// solved for c:
//   c ≡ store.ApplyInverse(store_elem, query.Apply(query_elem, q)).
// q is full-rank (independent constraints) and the group action is invertible,
// so the image is full-rank; we still drop any leading zero rows of the
// reversed RREF so the key matches the certificate's stored representatives.
template <class Problem>
Constraints<Problem::kP, Problem::kM, Problem::kNA> CanonicalFromWitness(
    const typename Problem::SymmetryGroup &group,
    const Constraints<Problem::kP, Problem::kM, Problem::kNA> &q,
    typename Problem::SymmetryGroup::QuerySet::Elem query_elem,
    typename Problem::SymmetryGroup::StoreSet::Elem store_elem) {
  constexpr int NA = Problem::kNA;
  constexpr int P = Problem::kP;
  constexpr int M = Problem::kM;
  using Vec = GFVec<P, M, NA>;
  Constraints<P, M, NA> image;
  image.reserve(q.size());
  for (const Vec v : q) {
    image.push_back(
        group.store.ApplyInverse(store_elem, group.query.Apply(query_elem, v)));
  }
  const int rank = GaussJordanRREF<P, M, NA>(&image);
  image.erase(image.begin(), image.end() - rank);
  return image;
}

template <class Problem> class RankLowerBoundBacktrackingVerifier {
public:
  static constexpr int NA = Problem::kNA;
  static constexpr int P = Problem::kP;
  static constexpr int M = Problem::kM;
  using Vec = GFVec<P, M, NA>;
  using SymmetryGroup = typename Problem::SymmetryGroup;
  using QueryElem = typename SymmetryGroup::QuerySet::Elem;
  using StoreElem = typename SymmetryGroup::StoreSet::Elem;

  static void Verify(const Constraints<P, M, NA> &constraints,
                     int rank_lower_bound, const BacktrackingProof &proof,
                     const SymmetryGroup &group, const RankMap<Problem> &map) {
    RankLowerBoundBacktrackingVerifier(constraints, rank_lower_bound, proof,
                                       group, map)
        .Verify();
  }

private:
  RankLowerBoundBacktrackingVerifier(const Constraints<P, M, NA> &constraints,
                                     int rank_lower_bound,
                                     const BacktrackingProof &proof,
                                     const SymmetryGroup &group,
                                     const RankMap<Problem> &map)
      : base_constraints_(constraints), rank_lower_bound_(rank_lower_bound),
        proof_(proof), group_(group), map_(map) {
    // Enumerate the same minimal constraints the prover does (see
    // rank_lower_bound_backtracking.h for the F_q minimality predicate).
    static_assert(NA < 32, "Backtracking DFS path mask is uint32_t");
    const uint64_t num_candidates = IntPow(P, M * NA);
    std::vector<int> base_pivots;
    base_pivots.reserve(base_constraints_.size());
    for (const auto &b : base_constraints_) {
      base_pivots.push_back(b.LeadingNonzeroIdx());
    }
    for (uint64_t value = 1; value < num_candidates; ++value) {
      const Vec constraint = DecodeGFVec<P, M, NA>(value);
      if (constraint.LeadingNonzero() != GF<P, M>::One()) {
        continue;
      }
      bool is_minimal = true;
      for (int p_idx : base_pivots) {
        if (constraint[p_idx] != GF<P, M>::Zero()) {
          is_minimal = false;
          break;
        }
      }
      if (is_minimal) {
        minimal_constraints_.push_back(constraint);
      }
    }
    CHECK(!minimal_constraints_.empty());
  }

  void Verify() const {
    size_t proof_index = 0;
    Constraints<P, M, NA> dfs_constraints;
    Verify(static_cast<int>(minimal_constraints_.size()) - 1, &dfs_constraints,
           &proof_index);
    CHECK_EQ(proof_index, proof_.Size());
  }

  void Verify(int max_constraint_idx, Constraints<P, M, NA> *dfs_constraints,
              size_t *proof_index) const {
    CHECK_LT(*proof_index, proof_.Size());
    const size_t proof_dfs_constraints_size =
        proof_.dfs_constraints_size_array[*proof_index];
    CHECK_LE(dfs_constraints->size(), proof_dfs_constraints_size);
    if (dfs_constraints->size() == proof_dfs_constraints_size) {
      const uint32_t mask = proof_.mask_array[*proof_index];
      const auto query_elem =
          static_cast<QueryElem>(proof_.query_elem_array[*proof_index]);
      const auto store_elem =
          static_cast<StoreElem>(proof_.store_elem_array[*proof_index]);
      Constraints<P, M, NA> extended = base_constraints_;
      for (int i = 0; i < static_cast<int>(dfs_constraints->size()); ++i) {
        if (mask & (uint32_t(1) << i)) {
          extended.push_back((*dfs_constraints)[i]);
        }
      }
      CHECK_GT(extended.size(), base_constraints_.size());
      const Constraints<P, M, NA> canonical = CanonicalFromWitness<Problem>(
          group_, extended, query_elem, store_elem);
      const auto it = map_.find(canonical);
      CHECK(it != map_.end())
          << "backtracking: canonical orbit not found in map";
      const int proved = std::popcount(mask) + it->second;
      CHECK_LE(rank_lower_bound_, proved);
      ++*proof_index;
      return;
    }

    for (int i = 0; i <= max_constraint_idx; ++i) {
      dfs_constraints->push_back(minimal_constraints_[i]);
      Verify(i, dfs_constraints, proof_index);
      dfs_constraints->pop_back();
    }
  }

  const Constraints<P, M, NA> &base_constraints_;
  int rank_lower_bound_ = 0;
  const BacktrackingProof &proof_;
  const SymmetryGroup &group_;
  const RankMap<Problem> &map_;
  std::vector<Vec> minimal_constraints_;
};
