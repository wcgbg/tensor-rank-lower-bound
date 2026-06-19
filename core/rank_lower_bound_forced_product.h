#pragma once

// Forced Product lower bound (Hopcroft–Kerr Lemma 2).
//
// Lemma 2: let F = {f_0, …, f_{n−1}} be expressions where f_0, …, f_{k−1} are
// independent and each is a single product. If F is computable with p
// multiplications, then there is an algorithm for F with p multiplications in
// which k of them are exactly f_0, …, f_{k−1}.
//
// Applied to the A-mode: each a-coordinate i contributes the bc-slice T[i] (an
// NB × NC matrix). The rank-1 slices are "single products" and, once a maximal
// independent set of r1 of them is fixed as forced products, the remaining
// rank must cover the higher-rank slices plus every 𝔽_q-combination of the
// forced products folded into them. We search all q^(r2p·r1) combinations and
// take the minimum flattening lower bound, then add back the r1 forced rows.
//
// SOUNDNESS NOTE: the lemma gives rank(T) = r1 + min over ALL combinations c of
// rank(T_c). We lower-bound each rank(T_c) by a flatten bound and take the min,
// so the min MUST range over every combination. The coefficients live in 𝔽_q
// (q = P^M), so for M ≥ 2 each coefficient ranges over all q field elements —
// enumerating only the prime subfield 𝔽_P would skip combinations and could
// report a min larger than the true rank (an unsound lower bound).
//
// Split out of the matrix-mult rank_lower_bound_basic_technics.h. Only the
// A-loop (RankLowerBoundForcedProductA) lives here; the three-position wrapper
// that records which projection won lives in rank_lower_bound_computer.h.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <vector>

#include <ng-log/logging.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include "core/dynamic_matrix.h"
#include "core/gf.h"
#include "core/rank_lower_bound_flatten.h"
#include "core/tensor.h"

// Returns a rank lower bound from the Forced Product technique on the A-mode,
// or 0 when the technique does not apply (no rank-1 slices, or the search space
// is too large to enumerate). Never returns more than max(known_lower_bound,…);
// the caller maxes the three cyclic positions.
//
// Each "coefficient position" can take any value in the full field 𝔽_q, so the
// enumeration covers q^bit_width combinations (q = P^M). The slice-update is
// `tensor_t[i][j][k] += digit · r1_bc[...]` over 𝔽_q, where `digit` ranges over
// every field element. The P=2, M=1 prime field is special-cased to a
// single-bit XOR fast path; every other (P, M) uses the general field path.
template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
int RankLowerBoundForcedProductA(const Tensor<P, M, NA, NB, NC> &tensor,
                                 int known_lower_bound = 0,
                                 int max_iterations_log2 = 64) {
  using Field = GF<P, M>;
  // Classify each a-slice by its NB × NC rank.
  std::vector<int> a_to_rank_bc(NA, 0);
  int r1_count = 0;
  for (std::size_t i = 0; i < NA; ++i) {
    DynamicMatrix<P, M> bc_matrix(tensor[i]);
    int rank = bc_matrix.Rank();
    a_to_rank_bc[i] = rank;
    if (rank == 1) {
      ++r1_count;
    }
  }
  if (r1_count == 0) {
    return 0;
  }

  // Partition slices: a maximal independent set of rank-1 slices becomes the
  // forced-product rows (r1_bc_collection); everything else (rank ≥ 2, plus the
  // rank-1 slices that were dependent on already-chosen ones) goes into r2p.
  Tensor<P, M, NA, NB, NC> r2p = {};
  int r2p_size0 = 0;
  DynamicMatrix<P, M> r1_bc_collection(0, NB * NC);
  int r1_bc_rows = 0;
  for (std::size_t i = 0; i < NA; ++i) {
    if (a_to_rank_bc[i] == 0) {
      // drop: zero slice contributes nothing
    } else if (a_to_rank_bc[i] == 1) {
      r1_bc_collection.ResizeRows(r1_bc_rows + 1);
      for (std::size_t j = 0; j < NB; ++j) {
        for (std::size_t k = 0; k < NC; ++k) {
          r1_bc_collection(r1_bc_rows, static_cast<int>(j * NC + k)) =
              tensor[i][j][k];
        }
      }
      if (r1_bc_collection.Rank() == r1_bc_rows + 1) {
        ++r1_bc_rows;
      } else {
        r1_bc_collection.ResizeRows(r1_bc_rows);
        r2p[r2p_size0] = tensor[i];
        ++r2p_size0;
      }
    } else {
      r2p[r2p_size0] = tensor[i];
      ++r2p_size0;
    }
  }

  constexpr int kQ = Field::kQ; // q = P^M; coefficients range over all of 𝔽_q.
  const int bit_width = r2p_size0 * r1_bc_rows;
  // num_iterations = kQ^bit_width. Cancel (return 0) when it exceeds the
  // enumeration budget. bit_width can be O(NA^2), so kQ^bit_width can be
  // astronomically larger than uint64 — we must NOT hand that to IntPow,
  // which now treats overflow as fatal.
  CHECK_GE(max_iterations_log2, 0);
  CHECK_LE(max_iterations_log2, 64);
  uint64_t num_iterations = 1;
  for (int i = 0; i < bit_width; ++i) {
    num_iterations *= static_cast<uint64_t>(kQ);
    if (max_iterations_log2 < 64 &&
        num_iterations > (uint64_t{1} << max_iterations_log2)) {
      LOG(WARNING) << "Cancel Forced Product. num_iterations=" << kQ << "^("
                   << r2p_size0 << "*" << r1_bc_rows << ") exceeds 2^"
                   << max_iterations_log2;
      return 0;
    }
  }
  if (num_iterations > (uint64_t{1} << 24)) {
    LOG(INFO) << "Forced Product. num_iterations=" << kQ << "^(" << r2p_size0
              << "*" << r1_bc_rows << ")=" << num_iterations;
  }

  const auto start_time = std::chrono::steady_clock::now();

  const Field *r1_bc_data = r1_bc_collection.data();
  CHECK_NOTNULL(r1_bc_data);
  // Scrambler: a large prime coprime to num_iterations gives a pseudo-random
  // visit order so an early break is likely to hit a low-rank witness fast.
  // For F_2 this matches the legacy code (kPrime modulo 2^32).
  constexpr uint64_t kPrime = 73074167;

  std::atomic<bool> early_break{false};
  std::atomic<uint64_t> progress{0};
  int rank_lower_bound = tbb::parallel_reduce(
      tbb::blocked_range<uint64_t>(0, num_iterations),
      std::numeric_limits<int>::max(),
      [&](const tbb::blocked_range<uint64_t> &range, int init) {
        int local_min = init;
        for (uint64_t t = range.begin(); t != range.end(); ++t) {
          uint64_t local_progress =
              progress.fetch_add(1, std::memory_order_relaxed);
          if (local_progress > 0 && local_progress % (uint64_t(1) << 24) == 0) {
            LOG(INFO) << std::format(
                "Progress: {}/{} = {:.2f}%", local_progress, num_iterations,
                static_cast<double>(local_progress) / num_iterations * 100.0);
          }
          if (early_break) {
            break;
          }
          // Scrambled visit order.
          const uint64_t scrambled = (t * kPrime) % num_iterations;
          Tensor<P, M, NA, NB, NC> tensor_t = r2p;
          if constexpr (P == 2 && M == 1) {
            // F_2 fast path: each coefficient is a bit, contribution is XOR.
            // Over GF<2,1>, `+=` is XOR.
            for (int bit_idx = 0; bit_idx < bit_width; ++bit_idx) {
              if (((scrambled >> bit_idx) & 1) == 0) {
                continue;
              }
              int r1_idx = bit_idx % r1_bc_rows; // forced-product row
              int i = bit_idx / r1_bc_rows;      // r2p slice
              for (std::size_t j = 0; j < NB; ++j) {
                for (std::size_t k = 0; k < NC; ++k) {
                  tensor_t[i][j][k] +=
                      r1_bc_data[r1_idx * (NB * NC) + j * NC + k];
                }
              }
            }
          } else {
            // General 𝔽_q path: read scrambled as a base-q numeral; each digit
            // is a field index in [0, q) (the GF<P,M> `value` IS the index), so
            // the coefficient ranges over every element of 𝔽_q. Contribution is
            // digit · r1_bc_data over 𝔽_q.
            uint64_t v = scrambled;
            for (int bit_idx = 0; bit_idx < bit_width; ++bit_idx) {
              const Field digit{static_cast<uint8_t>(v % kQ)};
              v /= kQ;
              if (digit == Field::Zero()) {
                continue;
              }
              int r1_idx = bit_idx % r1_bc_rows;
              int i = bit_idx / r1_bc_rows;
              for (std::size_t j = 0; j < NB; ++j) {
                for (std::size_t k = 0; k < NC; ++k) {
                  const Field r1_val =
                      r1_bc_data[r1_idx * (NB * NC) + j * NC + k];
                  tensor_t[i][j][k] = tensor_t[i][j][k] + digit * r1_val;
                }
              }
            }
          }
          int remaining = RankLowerBoundFlatten<P, M, NA, NB, NC>(
              tensor_t, local_min - r1_bc_rows);
          local_min = std::min(local_min, r1_bc_rows + remaining);
          if (local_min <= known_lower_bound) {
            early_break = true;
          }
        }
        return local_min;
      },
      [](int a, int b) { return std::min(a, b); });

  const auto end_time = std::chrono::steady_clock::now();
  const auto duration_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                            start_time)
          .count();
  if (duration_ms > 1000) {
    LOG(INFO) << std::format("RankLowerBoundForcedProductA. rank={} dur={:.2f}",
                             rank_lower_bound, duration_ms / 1000.0);
  }

  // An early break means we proved the minimum can't beat known_lower_bound, so
  // this position yields no improvement.
  if (early_break.load()) {
    return known_lower_bound;
  }
  return std::max(known_lower_bound, rank_lower_bound);
}
