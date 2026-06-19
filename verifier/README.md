# Verifier — the trusted checker

The correctness of every rank **lower** bound this project reports rests on this
directory. `verifier_main` reads a `pb::Certificate` (plus its sibling `.btp`
backtracking-proof archive) and independently re-checks every orbit's proof. If
it prints `OK`, the proven lower bound is sound; any failed check aborts via
`CHECK`. The verifier is deliberately far simpler than the search/prover that
*produced* the certificate, which is why it is the thing to audit.

## What you must trust (the trust base)

To trust a verified lower bound you only need to read:

- **`verifier/`** — this directory: the checker.
  - `verifier.h` (`core::VerifyRankLowerBound`) — the NA→0 dimension sweep and
    per-orbit proof dispatch.
  - `backtracking_verifier.h` — the plain `RankMap` (one canonical→rank entry
    per orbit) and the backtracking-trace replay.
  - `verifier_main.cc` — the entry point.
- **`core/`** — trusted primitives: field arithmetic (`gf`, `gf_vec`,
  `bit_vec`, `math_utils`, `dynamic_matrix`), `tensor`, `constraints`, the
  certificate proto + `proto_io`, the `backtracking_proof` codec, and the
  `rank_lower_bound_flatten` / `rank_lower_bound_forced_product` recomputers.
  `core/symmetry.h` is the `SymmetryGroupConcept` the problems implement and the
  verifier relies on (via witness reconstruction).
- **`problems/`** — the tensor and symmetry-group definitions
  (`Problem::MakeTensor`, `Problem::SymmetryGroup`). The verifier rebuilds the
  problem tensor from these, so they are part of the trust base.

## What you can ignore

`search/` (orbit enumeration, the meet-in-the-middle `OrbitMap`, the top-down DP
computer, and the backtracking-DFS / degenerate proof *generators*) and
`upper_bound/` (the flip-graph rank *upper* bound) are **not** part of the trust
base. They produce certificates; they do not vouch for them. A bug there can at
worst make the verifier reject — it cannot make it wrongly accept.

## What the verifier checks, per orbit

It repeats the prover's dimension sweep NA→0, but with a plain `RankMap` instead
of the prover's Store-image-expanded `OrbitMap`. For each orbit it re-checks the
recorded `RankLowerBoundProof`:

- **flatten** / **forced product** — by recomputation from the tensor.
- **degenerate reduction** — reconstructs the canonical form of the one-larger
  constraint set from the proof's `(query_elem, store_elem)` witness
  (`c ≡ store.ApplyInverse(store, query.Apply(query, q))`) and looks it up in the
  map; a wrong witness yields a missing key → `CHECK` failure.
- **backtracking** — replays the orbit's gzip-compressed blob from the once-loaded
  `.btp` archive via `backtracking_verifier.h`.

## Running it

A plain build uses the committed default problem in `chosen_problem.h`. To verify
a certificate for a specific problem, inject it with a build define (scoped to the
`*_main.cc` translation units):

The family is chosen with one of `CP_MATRIX` / `CP_CYCLIC` / `CP_FULL` /
`CP_EXTFIELD` and the numeric parameters as separate defines (kept comma-free so
they survive Bazel's comma-splitting of copt option lists):

```bash
bazel build --config=opt //verifier:verifier_main \
  --per_file_copt='.*_main\.cc@-DCP_MATRIX,-DCP_P=2,-DCP_M=1,-DCP_N0=3,-DCP_N1=3,-DCP_N2=3'
bazel-bin/verifier/verifier_main certs/matrix/cert_matrix_q02_n333.pb.txt
```

`run.py` passes these defines automatically (`run.py verify certs/matrix/cert_matrix_q02_n333.pb.txt`).
