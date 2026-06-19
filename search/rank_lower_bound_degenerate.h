#pragma once

// Degenerate Reduction lower bound.
//
// If appending one extra independent constraint to the current subspace lands
// in an orbit whose rank lower bound is r, then the current (less constrained)
// tensor also has rank ≥ r: constraining can only lower rank. So we sweep every
// possible extra constraint, keep those that stay linearly independent, look
// the resulting (one-larger) subspace up in the OrbitMap, and take the best.
//
// Faithful port of rank_search/rank_lower_bound_computer.h's
// RankLowerBoundDegenerate, with the matrix data type swapped for the NA-wide
// GFVec row and the (transpose, gl_left, gl_right) witness replaced by the
// OrbitMap's (query_elem, store_elem) witness.

#include <cstdint>
#include <utility>

#include <ng-log/logging.h>

#include "core/certificate.pb.h"
#include "core/constraints.h"
#include "core/gf.h"
#include "core/gf_vec.h"
#include "search/orbit_map.h"

template <class Problem>
std::pair<int, pb::DegenerateProof> RankLowerBoundDegenerate(
    const OrbitMap<Problem> &orbit_map,
    Constraints<Problem::kP, Problem::kM, Problem::kNA> constraints) {
  constexpr int NA = Problem::kNA;
  constexpr int P = Problem::kP;
  constexpr int M = Problem::kM;
  using Row = GFVec<P, M, NA>;
  using QueryElem = typename Problem::SymmetryGroup::QuerySet::Elem;
  using StoreElem = typename Problem::SymmetryGroup::StoreSet::Elem;

  int rank_lower_bound = 0;
  pb::DegenerateProof degenerate_proof;
  // Enumerate every q^NA = P^(M*NA) candidate row, decoded from an integer
  // counter. For (P,M)≠(2,1) we additionally require leading-1 (the canonical
  // RREF normalization); for (P,M)=(2,1) every nonzero row already has
  // leading 1.
  const uint64_t num_candidates = IntPow(P, M * NA);
  for (uint64_t r = 1; r < num_candidates; ++r) {
    Row constraint = DecodeGFVec<P, M, NA>(r);
    if (constraint.LeadingNonzero() != GF<P, M>::One()) {
      continue;
    }
    constraints.push_back(constraint);
    if (IsLinearIndependentRREF<P, M, NA>(constraints)) {
      QueryElem query_elem{};
      StoreElem store_elem{};
      const int found = orbit_map.Get(constraints, &query_elem, &store_elem);
      if (found > rank_lower_bound) {
        rank_lower_bound = found;
        const uint64_t encoded = EncodeGFVec<P, M, NA>(constraint);
        CHECK_LT(encoded, uint64_t{1} << 32)
            << "extra_constraint overflows fixed32 (P=" << P << ", M=" << M
            << ", NA=" << NA << "); widen the proto field if you hit this.";
        degenerate_proof.set_extra_constraint(static_cast<uint32_t>(encoded));
        // The witness elements double as their group-enumeration indices (the
        // opaque fixed32 the proto stores); this holds for groups whose
        // query.At/store.At are the identity on indices, e.g. cyclic.
        degenerate_proof.mutable_transformation()->set_query_elem(
            static_cast<uint32_t>(query_elem));
        degenerate_proof.mutable_transformation()->set_store_elem(
            static_cast<uint32_t>(store_elem));
      }
    }
    constraints.pop_back();
  }
  return {rank_lower_bound, degenerate_proof};
}
