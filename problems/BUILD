load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

# Compile-time problem selection shared by the prover binaries (search/,
# upper_bound/) and the verifier (verifier/). A selector over the problem
# definitions in this package's subdirectories; the problem is chosen via
# -DCP_* build defines — see chosen_problem.h.
cc_library(
    name = "chosen_problem",
    hdrs = ["chosen_problem.h"],
    deps = [
        "//problems/cyclic:problem",
        "//problems/extfield:problem",
        "//problems/full:problem",
        "//problems/matrix:problem",
    ],
)
