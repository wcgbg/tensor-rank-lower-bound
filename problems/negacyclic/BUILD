load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "tensor",
    hdrs = ["tensor.h"],
    deps = [
        "//core:gf",
        "//core:tensor",
    ],
)

cc_test(
    name = "tensor_test",
    srcs = ["tensor_test.cc"],
    deps = [
        ":tensor",
        "//core:gf",
        "//core:tensor",
        "//core:tensor_utils",
        "//problems/cyclic:tensor",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "symmetry",
    hdrs = ["symmetry.h"],
    deps = [
        "//problems/poly_quotient:symmetry",
    ],
)

cc_library(
    name = "problem",
    hdrs = ["problem.h"],
    deps = [
        ":symmetry",
        ":tensor",
        "//core:bit_vec",
        "//core:constraints",
        "//core:gf_vec",
        "//core:tensor",
    ],
)

cc_test(
    name = "symmetry_test",
    srcs = ["symmetry_test.cc"],
    deps = [
        ":symmetry",
        ":tensor",
        "//core:bit_vec",
        "//core:constraints",
        "//core:gf",
        "//core:gf_vec",
        "//core:rank_lower_bound_flatten",
        "//core:tensor",
        "//core:tensor_isomorphism",
        "//core:tensor_utils",
        "//problems/cyclic:symmetry",
        "@googletest//:gtest_main",
    ],
)
