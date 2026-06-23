load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "irreducibles",
    hdrs = ["irreducibles.h"],
    deps = [
        "//core:gf",
        "//core:gf_vec",
    ],
)

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
        ":irreducibles",
        ":tensor",
        "//core:gf",
        "//core:gf_vec",
        "//core:tensor",
        "//core:tensor_utils",
        "//problems/cyclic:tensor",
        "//problems/truncated:tensor",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "symmetry",
    hdrs = ["symmetry.h"],
    deps = [
        ":irreducibles",
        "//core:bit_vec",
        "//core:gf",
        "//core:gf_vec",
        "//core:constraints",
    ],
)

cc_test(
    name = "symmetry_test",
    srcs = ["symmetry_test.cc"],
    deps = [
        ":irreducibles",
        ":symmetry",
        ":tensor",
        "//core:bit_vec",
        "//core:constraints",
        "//core:tensor",
        "//core:tensor_isomorphism",
        "//core:tensor_utils",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "problem",
    hdrs = ["problem.h"],
    deps = [
        ":irreducibles",
        ":symmetry",
        ":tensor",
        "//core:bit_vec",
        "//core:gf_vec",
        "//core:constraints",
        "//core:tensor",
    ],
)
