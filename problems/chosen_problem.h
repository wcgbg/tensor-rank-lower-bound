#pragma once

// Compile-time problem selection, shared by the prover binaries (in search/ and
// upper_bound/) and the verifier (in verifier/). It lives in problems/ — the
// common dependency of all of them — as a selector over the problem definitions
// rather than something tied to either the prover or the verifier side.
//
// The problem is chosen at *build time* via preprocessor defines, NOT by
// editing this file. Select a family with one of CP_MATRIX / CP_CYCLIC /
// CP_FULL / CP_EXTFIELD and pass the numeric parameters as separate defines.
// They are kept comma-free on purpose: Bazel splits copt values on commas, so
// the commas below separate distinct -D options (not template arguments).
// Examples:
//
//   # matrix <2,1,3,3,3>
//   bazel build //verifier:verifier_main \
//     --per_file_copt='.*_main\.cc@-DCP_MATRIX,-DCP_P=2,-DCP_M=1,-DCP_N0=3,-DCP_N1=3,-DCP_N2=3'
//   # cyclic <2,1,7>
//   bazel build //search:rank_lower_bound_main \
//     --per_file_copt='.*_main\.cc@-DCP_CYCLIC,-DCP_P=2,-DCP_M=1,-DCP_N=7'
//
// chosen_problem.h is included only by the *_main.cc binaries, so scoping the
// define to those translation units avoids invalidating the rest of the build
// cache (core/, problems/, external deps) when the problem changes. run.py
// passes these flags automatically. A plain `bazel build //...` with no
// override uses the committed default below, so `//...` builds and tests keep
// working.
//
// The polynomial-multiplication problems are parameterized as `<P, M, N>`:
//   - P: characteristic of the base field (P prime).
//   - M: extension degree of the base field. The field is 𝔽_q with q = P^M.
//   - N: polynomial dimension.
// The matrix-multiplication family takes three size parameters instead of one:
// `matrix::Problem<P, M, N0, N1, N2>` (the ⟨N0, N1, N2⟩ format).
//
// rank_upper_bound_main is hard-coded to (P, M) == (2, 1) (the F_2 flip-graph
// search). For q ≠ 2 runs skip that stage entirely — see run.py.

#include "problems/cyclic/problem.h"
#include "problems/extfield/problem.h"
#include "problems/full/problem.h"
#include "problems/matrix/problem.h"

#if defined(CP_MATRIX)
#define CHOSEN_PROBLEM matrix::Problem<CP_P, CP_M, CP_N0, CP_N1, CP_N2>
#elif defined(CP_CYCLIC)
#define CHOSEN_PROBLEM cyclic::Problem<CP_P, CP_M, CP_N>
#elif defined(CP_FULL)
#define CHOSEN_PROBLEM full::Problem<CP_P, CP_M, CP_N>
#elif defined(CP_EXTFIELD)
#define CHOSEN_PROBLEM extfield::Problem<CP_P, CP_M, CP_N>
#else
// Committed default, used when no CP_* define is passed (e.g. plain
// `bazel build //...`). Matches matrix::Problem<2, 1, 3, 3, 3>.
#define CHOSEN_PROBLEM matrix::Problem<2, 1, 3, 3, 3>
#endif

// Declared at global scope (no namespace): it is a single project-wide
// build-time selection, used by the *_main.cc binaries as `::ChosenProblem`.
using ChosenProblem = CHOSEN_PROBLEM;
