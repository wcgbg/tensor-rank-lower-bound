#pragma once

// Gauss-Jordan elimination over 𝔽_q with q = P^M.
//
// Two entry points live here:
//   - GaussJordanEliminationF2<N> for the (P, M) == (2, 1) bit-packed hot path
//     (operates on the BitVec<N> inside each GFVec<2, 1, N>);
//   - GaussJordanEliminationFq<P, M, N> for every other (P, M), operating on
//     the packed-digit GFVec<P, M, N> rows.
//
// Both produce full reduced row-echelon form (RREF):
//   - column order reversed (highest-coordinate pivot first);
//   - row order reversed (zero rows park at the front), so erasing the
//     leading-zero rows gives the canonical key the OrbitMap stores;
//   - each pivot row is normalized to leading 1 — for F_q this is a scalar
//     multiplication by the pivot's inverse; for F₂ the pivot is already 1;
//   - entries above and below each pivot are eliminated (full RREF, not just
//     row echelon).
//
// Returns the number of nonzero rows (the rank). core/constraints.h provides
// a dispatching wrapper (GaussJordanRREF) that picks the right entry point.

#include <vector>

#include "core/bit_vec.h"
#include "core/gf.h"
#include "core/gf_vec.h"

// F₂ Gauss-Jordan elimination on the bit-packed (P, M) == (2, 1) row
// representation. Accesses each row's underlying BitVec via `.data` so the
// algorithm operates directly on the GFVec<2, 1, N> wrapper, preserving the
// XOR/popcount hot path.
//
// Outputs are full RREF in the column-reversed convention:
//   - row order reversed (zero rows park at the front),
//   - each pivot row is leading-1 (already true in F₂),
//   - above-and-below elimination (full Jordan).
// Returns the number of nonzero rows (the rank).
template <int N>
int GaussJordanEliminationF2(std::vector<GFVec<2, 1, N>> *matrix) {
  using BV = BitVec<N>;
  const int size = static_cast<int>(matrix->size());
  if (size == 0 || N <= 0) {
    return 0;
  }
  GFVec<2, 1, N> *const data = matrix->data();

  int rank = 0;
  int current_row = size - 1;
  int pivot_col = N - 1;

  while (current_row >= 0 && pivot_col >= 0) {
    const BV col_mask = static_cast<BV>(BV{1} << pivot_col);
    int pivot_row = -1;
    for (int r = current_row; r >= 0; --r) {
      if (data[r].data & col_mask) {
        pivot_row = r;
        break;
      }
    }
    if (pivot_row == -1) {
      --pivot_col;
      continue;
    }

    if (pivot_row != current_row) {
      std::swap(data[current_row], data[pivot_row]);
    }

    for (int r = 0; r < size; ++r) {
      if (r != current_row && (data[r].data & col_mask)) {
        data[r].data = static_cast<BV>(data[r].data ^ data[current_row].data);
      }
    }

    ++rank;
    --current_row;
    --pivot_col;
  }

  return rank;
}

template <int P, int M, int N>
int GaussJordanEliminationFq(std::vector<GFVec<P, M, N>> *matrix) {
  static_assert(!(P == 2 && M == 1),
                "GaussJordanEliminationFq is for (P,M) != (2,1); use "
                "GaussJordanEliminationF2 for the F₂ BitVec hot path");
  using Vec = GFVec<P, M, N>;
  using Field = GF<P, M>;

  const int size = static_cast<int>(matrix->size());
  if (size == 0 || N <= 0) {
    return 0;
  }
  Vec *const data = matrix->data();

  int rank = 0;
  int current_row = size - 1;
  int pivot_col = N - 1;

  const Field zero = Field::Zero();
  const Field one = Field::One();

  while (current_row >= 0 && pivot_col >= 0) {
    // Find a pivot row: any row in [0, current_row] with a nonzero entry in
    // the pivot column.
    int pivot_row = -1;
    for (int r = current_row; r >= 0; --r) {
      if (data[r][pivot_col] != zero) {
        pivot_row = r;
        break;
      }
    }
    if (pivot_row == -1) {
      --pivot_col;
      continue;
    }

    if (pivot_row != current_row) {
      std::swap(data[current_row], data[pivot_row]);
    }

    // Normalize the pivot row to leading 1 by scaling by the inverse of the
    // pivot. Cheap: a single field-scalar multiplication.
    const Field pivot = data[current_row][pivot_col];
    if (pivot != one) {
      data[current_row] = pivot.Inverse() * data[current_row];
    }

    // Eliminate above and below the pivot: for each other row r with a
    // nonzero entry c in this column, do row[r] -= c * row[current_row].
    for (int r = 0; r < size; ++r) {
      if (r == current_row) {
        continue;
      }
      const Field c = data[r][pivot_col];
      if (c == zero) {
        continue;
      }
      data[r] = data[r] - (c * data[current_row]);
    }

    ++rank;
    --current_row;
    --pivot_col;
  }

  return rank;
}

template <int P, int M, int N>
bool IsLinearIndependentFq(std::vector<GFVec<P, M, N>> matrix) {
  return GaussJordanEliminationFq<P, M, N>(&matrix) ==
         static_cast<int>(matrix.size());
}
