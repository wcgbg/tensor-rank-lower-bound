load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

# Shared A-side symmetry-group engine for the three polynomial-quotient problem
# families that are multiplication in F_q[x]/(x^N - gamma): cyclic (gamma=1),
# truncated (gamma=0), and negacyclic (gamma=-1). Parameterised by the reduction
# constant via the kGamma template argument; each family aliases it in its own
# symmetry.h.
cc_library(
    name = "symmetry",
    hdrs = ["symmetry.h"],
    deps = [
        "//core:bit_vec",
        "//core:gf",
        "//core:gf_vec",
        "@boost.unordered",
        "@ng-log//:ng-log",
        "@onetbb//:tbb",
    ],
)
