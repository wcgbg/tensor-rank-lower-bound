# The rank-lower-bound verifier — the trusted checker. Together with //core
# (trusted primitives) and //problems (tensor + symmetry definitions, including
# //problems:chosen_problem), this is the complete trust base for the project's rank
# lower bounds. See verifier/README.md. It does NOT depend on //search or
# //upper_bound.
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "backtracking_verifier",
    hdrs = ["backtracking_verifier.h"],
    deps = [
        "//core:backtracking_proof",
        "//core:bit_vec",
        "//core:constraints",
        "@boost.unordered",
        "@ng-log//:ng-log",
    ],
)

cc_library(
    name = "verifier",
    hdrs = ["verifier.h"],
    deps = [
        ":backtracking_verifier",
        "//core:backtracking_proof",
        "//core:bit_vec",
        "//core:certificate_cc_proto",
        "//core:constraints",
        "//core:rank_lower_bound_flatten",
        "//core:rank_lower_bound_forced_product",
        "//core:tensor",
        "@boost.unordered",
        "@ng-log//:ng-log",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "verifier_main",
    srcs = ["verifier_main.cc"],
    malloc = "@mimalloc",
    deps = [
        ":verifier",
        "//core:backtracking_proof",
        "//core:certificate_cc_proto",
        "//core:proto_io",
        "//problems:chosen_problem",
        "@gflags",
        "@ng-log//:ng-log",
    ],
)

# Generates certificates at runtime (via //search) and verifies them, so the
# test — but not the trusted verifier library above — depends on //search.
cc_test(
    name = "verifier_test",
    srcs = ["verifier_test.cc"],
    deps = [
        ":verifier",
        "//core:backtracking_proof",
        "//core:certificate_cc_proto",
        "//problems/cyclic:problem",
        "//problems/extfield:problem",
        "//problems/full:problem",
        "//search:orbit_enumerator_slow",
        "//search:rank_lower_bound_computer",
        "@googletest//:gtest_main",
        "@onetbb//:tbb",
    ],
)
