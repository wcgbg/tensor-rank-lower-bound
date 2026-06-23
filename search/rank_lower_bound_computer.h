#pragma once

// Top-level dynamic-programming driver.
//
// Reads a Certificate (a list of orbit representatives with their current rank
// lower bounds) and iteratively improves each orbit's lower bound using the
// four techniques, sweeping subspace dimensions from most-constrained (NA, the
// full dual space) down to least-constrained (0, the empty subspace). Because a
// dimension's degenerate / backtracking lookups only ever consult orbits of
// strictly larger dimension — already processed and inserted into the OrbitMap
// — one top-down pass per dimension settles every orbit.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <ng-log/logging.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "core/backtracking_proof.h"
#include "core/bit_vec.h"
#include "core/certificate.pb.h"
#include "core/constraints.h"
#include "core/proto_io.h"
#include "core/rank_lower_bound_flatten.h"
#include "core/rank_lower_bound_forced_product.h"
#include "core/tensor.h"
#include "search/orbit_map.h"
#include "search/rank_lower_bound_backtracking.h"
#include "search/rank_lower_bound_degenerate.h"

struct ProcessOptions {
  bool basic_method = true; // Flatten + Forced Product
  bool degenerate_method = true;
  uint64_t backtracking_step_limit = std::numeric_limits<uint64_t>::max();
  size_t backtracking_max_map_size = 10'000'000;
  int dim_min = 0;
  int dim_max = std::numeric_limits<int>::max();
  int forced_product_max_iterations_log2 = 24;
};

// Three-position Forced Product: run RankLowerBoundForcedProductA on the tensor
// and its two cyclic transposes, recording which projection won. Mirrors
// rank_search/rank_lower_bound_computer.h's RankLowerBoundForcedProduct.
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
std::pair<int, pb::ForcedProductProof>
RankLowerBoundForcedProduct(const Tensor<P, M, NA, NB, NC> &tensor,
                            int known_lower_bound, int max_iterations_log2) {
  int rank0 = RankLowerBoundForcedProductA<P, M, NA, NB, NC>(
      tensor, known_lower_bound, max_iterations_log2);
  known_lower_bound = std::max(known_lower_bound, rank0);
  Tensor<P, M, NB, NC, NA> tensor1 = CyclicTranspose<P, M, NA, NB, NC>(tensor);
  int rank1 = RankLowerBoundForcedProductA<P, M, NB, NC, NA>(
      tensor1, known_lower_bound, max_iterations_log2);
  known_lower_bound = std::max(known_lower_bound, rank1);
  Tensor<P, M, NC, NA, NB> tensor2 = CyclicTranspose<P, M, NB, NC, NA>(tensor1);
  int rank2 = RankLowerBoundForcedProductA<P, M, NC, NA, NB>(
      tensor2, known_lower_bound, max_iterations_log2);
  known_lower_bound = std::max(known_lower_bound, rank2);

  pb::ForcedProductProof forced_product_proof;
  if (rank0 == known_lower_bound) {
    forced_product_proof.set_projection_type(0);
  } else if (rank1 == known_lower_bound) {
    forced_product_proof.set_projection_type(1);
  } else if (rank2 == known_lower_bound) {
    forced_product_proof.set_projection_type(2);
  } else {
    return {0, forced_product_proof};
  }
  return {known_lower_bound, forced_product_proof};
}

// Number of constraints (subspace dimension) encoded in a ConstrainedTensor.
template <class Problem> int NumConstraints(const pb::ConstrainedTensor &rt) {
  return static_cast<int>(
      rt.constraints().size() /
      sizeof(GFVec<Problem::kP, Problem::kM, Problem::kNA>));
}

// Process a single constrained-tensor orbit, returning (new_lb, proof, trace).
// `trace` is the winning backtracking DFS trace when a backtracking proof won
// (the caller stores it in the per-certificate archive), and empty otherwise.
template <class Problem>
std::tuple<int, pb::RankLowerBoundProof, BacktrackingProof>
ProcessOrbit(const pb::ConstrainedTensor &rt,
             const OrbitMap<Problem> &orbit_map,
             const ProcessOptions &options) {
  constexpr int NA = Problem::kNA;
  constexpr int NB = Problem::kNB;
  constexpr int NC = Problem::kNC;
  constexpr int P = Problem::kP;
  constexpr int M = Problem::kM;

  int rank_lower_bound = rt.has_rank_lower_bound() ? rt.rank_lower_bound() : -1;
  pb::RankLowerBoundProof proof;

  Constraints<P, M, NA> constraints =
      ConstraintsFromBytes<P, M, NA>(rt.constraints());
  Tensor<P, M, NA, NB, NC> tensor = ApplyConstraintsToTensor<P, M, NA, NB, NC>(
      constraints, Problem::MakeTensor());

  if (options.basic_method && !rt.has_rank_lower_bound()) {
    // Flatten.
    int flatten_rank = RankLowerBoundFlatten<P, M, NA, NB, NC>(tensor);
    if (flatten_rank > rank_lower_bound) {
      rank_lower_bound = flatten_rank;
      *proof.mutable_flatten_matrix_proof() = {};
    }
  }

  if (options.degenerate_method) {
    auto [degenerate_rank, degenerate_proof] =
        RankLowerBoundDegenerate<Problem>(orbit_map, constraints);
    if (degenerate_rank > rank_lower_bound) {
      rank_lower_bound = degenerate_rank;
      *proof.mutable_degenerate_proof() = std::move(degenerate_proof);
    }
  }

  if (options.basic_method && !rt.has_rank_lower_bound()) {
    // Forced Product.
    auto [forced_product_rank, forced_product_proof] =
        RankLowerBoundForcedProduct<P, M, NA, NB, NC>(
            tensor, rank_lower_bound,
            options.forced_product_max_iterations_log2);
    if (forced_product_rank > rank_lower_bound) {
      rank_lower_bound = forced_product_rank;
      *proof.mutable_forced_product_proof() = std::move(forced_product_proof);
    }
  }

  // `blob` ends up holding the last winning trace iff backtracking ultimately
  // won (it is only assigned alongside proof.backtracking_proof), so it stays
  // empty when no backtracking proof is recorded.
  BacktrackingProof blob;
  if (options.backtracking_step_limit > 0) {
    while (true) {
      auto [backtracking_rank, backtracking_proof, backtracking_blob] =
          RankLowerBoundBacktracking<Problem>::Search(
              constraints, orbit_map, rank_lower_bound,
              options.backtracking_step_limit,
              options.backtracking_max_map_size);
      if (backtracking_rank <= rank_lower_bound) {
        break;
      }
      rank_lower_bound = backtracking_rank;
      *proof.mutable_backtracking_proof() = std::move(backtracking_proof);
      blob = std::move(backtracking_blob);
    }
  }

  return {rank_lower_bound, std::move(proof), std::move(blob)};
}

// Process all orbits at a single subspace dimension, then insert their
// (updated) bounds into the OrbitMap so the next-smaller dimension can consult
// them. Returns true iff any orbit's lower bound improved. Mirrors
// ProcessOneRankLowerBound in the matrix-mult code.
template <class Problem>
bool ProcessOrbitsAtDim(int dim, const ProcessOptions &options,
                        BacktrackingProofArchive *archive,
                        pb::Certificate *certificate,
                        OrbitMap<Problem> *orbit_map, std::mt19937_64 *rng) {
  constexpr int NA = Problem::kNA;
  constexpr int P = Problem::kP;
  constexpr int M = Problem::kM;

  const auto iteration_start = std::chrono::steady_clock::now();

  std::vector<pb::ConstrainedTensor *> rts;
  for (int i = 0; i < certificate->constrained_tensors_size(); ++i) {
    pb::ConstrainedTensor *rt = certificate->mutable_constrained_tensors(i);
    if (NumConstraints<Problem>(*rt) == dim) {
      rts.push_back(rt);
    }
  }
  std::shuffle(rts.begin(), rts.end(), *rng);
  LOG(INFO) << "Processing dim=" << dim << ", count=" << rts.size();

  // Phase 1: compute the new bound for every orbit (read-only on the map). Each
  // task writes only its own disjoint slot, so the backtracking blob rides
  // along with no shared mutation; the archive itself is touched only in serial
  // phase 3.
  std::vector<std::tuple<int, pb::RankLowerBoundProof, BacktrackingProof>>
      rank_and_proof_list(rts.size());
  std::atomic<int> progress = 0;
  tbb::parallel_for(0, static_cast<int>(rts.size()), [&](int idx) {
    auto [rank, rank_proof, blob] =
        ProcessOrbit<Problem>(*rts[idx], *orbit_map, options);
    if (rank_proof.proof_case() != pb::RankLowerBoundProof::PROOF_NOT_SET) {
      LOG(INFO) << "Better LB for index=" << rts[idx]->index() << ": "
                << rts[idx]->rank_lower_bound() << "->" << rank;
    }
    rank_and_proof_list[idx] = {rank, std::move(rank_proof), std::move(blob)};

    int local_progress = progress.fetch_add(1, std::memory_order_relaxed) + 1;
    double progress_percentage = local_progress * 100.0 / rts.size();
    LOG_EVERY_T(INFO, 10) << std::format(
        "ProcessOrbitsAtDim P1: {} / {} = {:.2f}%", local_progress, rts.size(),
        progress_percentage);
  });

  LOG(INFO) << "ProcessOrbitsAtDim P2...";
  // Phase 2: insert every orbit's (final) bound into the map so smaller
  // dimensions can look it up. Done after phase 1 so no Set races a Get.
  tbb::parallel_for(
      tbb::blocked_range<int>(0, static_cast<int>(rts.size())),
      [&](const tbb::blocked_range<int> &range) {
        for (int i = range.begin(); i != range.end(); ++i) {
          orbit_map->Set(ConstraintsFromBytes<P, M, NA>(rts[i]->constraints()),
                         std::get<0>(rank_and_proof_list[i]));
        }
      });

  LOG(INFO) << "ProcessOrbitsAtDim P3...";
  // Phase 3: write the improved bounds and proofs back into the certificate,
  // and update the in-memory archive (serial; the archive is never touched in a
  // parallel_for).
  bool has_update = false;
  for (int i = 0; i < static_cast<int>(rts.size()); ++i) {
    auto &[new_rank, new_proof, new_blob] = rank_and_proof_list[i];
    if (new_proof.proof_case() == pb::RankLowerBoundProof::PROOF_NOT_SET) {
      continue;
    }
    pb::ConstrainedTensor *rt = rts[i];
    CHECK_LE(rt->rank_lower_bound(), new_rank);
    has_update = true;
    if (new_proof.has_backtracking_proof()) {
      archive->Set(static_cast<int>(rt->index()), std::move(new_blob));
    } else {
      // A non-backtracking proof won; drop any stale blob for this orbit.
      archive->Clear(static_cast<int>(rt->index()));
    }
    rt->set_rank_lower_bound(new_rank);
    *rt->mutable_rank_lower_bound_proof() = std::move(new_proof);
    if (rt->has_rank_upper_bound()) {
      CHECK_LE(rt->rank_lower_bound(), rt->rank_upper_bound()) << rt->index();
    }
  }

  const auto iteration_end = std::chrono::steady_clock::now();
  LOG(INFO)
      << "ProcessOrbitsAtDim done. duration=" << std::fixed
      << std::setprecision(1)
      << std::chrono::duration<double>(iteration_end - iteration_start).count();
  return has_update;
}

// Run the full DP: pre-seed the OrbitMap with any bounds already in the
// certificate, then sweep dimensions from NA down to 0, processing each and
// checkpointing the certificate to `output_path` after every dimension.
template <class Problem>
void ProcessOrbits(const ProcessOptions &options,
                   const std::string &output_path,
                   pb::Certificate *certificate) {
  constexpr int NA = Problem::kNA;
  constexpr int P = Problem::kP;
  constexpr int M = Problem::kM;

  const typename Problem::SymmetryGroup group;
  OrbitMap<Problem> orbit_map(&group);
  std::mt19937_64 rng;

  // Pre-seed: orbits with a known bound (e.g. from a prior run, or dimensions
  // excluded by dim_min/dim_max) must be in the map for the boundary
  // dimension's degenerate / backtracking lookups to find them.
  tbb::parallel_for(
      tbb::blocked_range<int>(0, certificate->constrained_tensors_size()),
      [&](const tbb::blocked_range<int> &range) {
        for (int i = range.begin(); i != range.end(); ++i) {
          const pb::ConstrainedTensor &rt = certificate->constrained_tensors(i);
          if (rt.has_rank_lower_bound()) {
            orbit_map.Set(ConstraintsFromBytes<P, M, NA>(rt.constraints()),
                          rt.rank_lower_bound());
          }
        }
      });

  // Backtracking traces live in one archive next to the certificate output.
  // With no output path there is nowhere to anchor it, so it is skipped. For
  // partial re-runs (dim_min/dim_max), load any existing archive so orbits
  // whose backtracking proofs are not reprocessed keep their traces.
  const std::string archive_path =
      output_path.empty() ? std::string{}
                          : GetBacktrackingProofArchivePath(output_path);
  BacktrackingProofArchive archive;
  if (!archive_path.empty() && std::filesystem::exists(archive_path)) {
    archive = BacktrackingProofArchive::Load(archive_path);
    CHECK_EQ(archive.Size(),
             static_cast<size_t>(certificate->constrained_tensors_size()));
  } else {
    archive.Resize(certificate->constrained_tensors_size());
  }

  for (int dim = NA; dim >= 0; --dim) {
    if (dim < options.dim_min || dim > options.dim_max) {
      continue;
    }
    ProcessOrbitsAtDim<Problem>(dim, options, &archive, certificate, &orbit_map,
                                &rng);
    if (!output_path.empty()) {
      WriteProtoToFile(*certificate, output_path);
      archive.Save(archive_path);
    }
  }
}
