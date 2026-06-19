load("@protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

# This package holds the trusted primitives shared by the verifier (//verifier)
# and the search/prover (//search): field arithmetic, vectors, tensors,
# constraints, the certificate proto + IO, the backtracking-proof codec, and the
# flatten / forced-product rank recomputers. It is part of the verifier trust
# base (see verifier/README.md). The complex search-only code lives in //search.
#
# symmetry.h (the SymmetryGroupConcept) stays here even though only //search
# includes it directly: it is the trusted interface contract the //problems
# symmetry groups implement and the verifier relies on via witness
# reconstruction (store.ApplyInverse).

proto_library(
    name = "certificate_proto",
    srcs = ["certificate.proto"],
)

cc_proto_library(
    name = "certificate_cc_proto",
    deps = [":certificate_proto"],
)

cc_library(
    name = "bit_vec",
    hdrs = ["bit_vec.h"],
)

cc_library(
    name = "gf",
    hdrs = ["gf.h"],
    deps = ["@ng-log//:ng-log"],
)

cc_test(
    name = "gf_test",
    srcs = ["gf_test.cc"],
    deps = [
        ":gf",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "gf_vec",
    hdrs = ["gf_vec.h"],
    deps = [
        ":bit_vec",
        ":gf",
    ],
)

cc_test(
    name = "gf_vec_test",
    srcs = ["gf_vec_test.cc"],
    deps = [
        ":gf_vec",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "math_utils",
    hdrs = ["math_utils.h"],
    deps = [
        ":bit_vec",
        ":gf",
        ":gf_vec",
        "@ng-log//:ng-log",
    ],
)

cc_test(
    name = "math_utils_test",
    srcs = ["math_utils_test.cc"],
    deps = [
        ":gf_vec",
        ":math_utils",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "tensor",
    hdrs = ["tensor.h"],
    deps = [
        ":gf",
    ],
)

cc_library(
    name = "tensor_utils",
    hdrs = ["tensor_utils.h"],
    deps = [
        ":gf",
        ":tensor",
    ],
)

cc_test(
    name = "tensor_utils_test",
    srcs = ["tensor_utils_test.cc"],
    deps = [
        ":tensor_utils",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "tensor_isomorphism",
    testonly = True,
    hdrs = ["tensor_isomorphism.h"],
    deps = [
        ":bit_vec",
        ":constraints",
        ":gf",
        ":gf_vec",
        ":tensor",
    ],
)

cc_test(
    name = "tensor_isomorphism_test",
    srcs = ["tensor_isomorphism_test.cc"],
    deps = [
        ":tensor_isomorphism",
        ":tensor_utils",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "constraints",
    hdrs = ["constraints.h"],
    deps = [
        ":bit_vec",
        ":gf_vec",
        ":math_utils",
        ":tensor",
        "@ng-log//:ng-log",
    ],
)

cc_test(
    name = "constraints_test",
    srcs = ["constraints_test.cc"],
    deps = [
        ":bit_vec",
        ":constraints",
        ":gf_vec",
        ":tensor",
        ":tensor_utils",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "symmetry",
    hdrs = ["symmetry.h"],
    deps = [
        ":bit_vec",
    ],
)

cc_library(
    name = "proto_io",
    hdrs = ["proto_io.h"],
    deps = [
        "@boost.iostreams",
        "@ng-log//:ng-log",
        "@protobuf",
    ],
)

cc_test(
    name = "proto_io_test",
    srcs = ["proto_io_test.cc"],
    # The golden lives with the rest of the test inputs in //search; see the
    # filegroup comment there.
    data = ["//search:testdata_proto_io_golden"],
    deps = [
        ":certificate_cc_proto",
        ":proto_io",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "rank_lower_bound_flatten",
    hdrs = ["rank_lower_bound_flatten.h"],
    deps = [
        ":dynamic_matrix",
        ":tensor",
    ],
)

cc_library(
    name = "rank_lower_bound_forced_product",
    hdrs = ["rank_lower_bound_forced_product.h"],
    deps = [
        ":dynamic_matrix",
        ":gf",
        ":rank_lower_bound_flatten",
        ":tensor",
        "@ng-log//:ng-log",
        "@onetbb//:tbb",
    ],
)

cc_library(
    name = "backtracking_proof",
    srcs = ["backtracking_proof.cc"],
    hdrs = ["backtracking_proof.h"],
    deps = [
        "@boost.iostreams",
        "@ng-log//:ng-log",
    ],
)

cc_test(
    name = "backtracking_proof_test",
    srcs = ["backtracking_proof_test.cc"],
    deps = [
        ":backtracking_proof",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "dynamic_matrix",
    srcs = ["dynamic_matrix.cc"],
    hdrs = ["dynamic_matrix.h"],
    deps = [
        ":gf",
        "@ng-log//:ng-log",
    ],
)

cc_test(
    name = "dynamic_matrix_test",
    srcs = ["dynamic_matrix_test.cc"],
    deps = [
        ":dynamic_matrix",
        "@googletest//:gtest_main",
    ],
)
