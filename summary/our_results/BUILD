load("@rules_cc//cc:defs.bzl", "cc_binary")

package(default_visibility = ["//visibility:public"])

# Tooling that summarizes cached certificates into the tables under
# summary/our_results.

cc_binary(
    name = "cert_summary_main",
    srcs = ["cert_summary_main.cc"],
    malloc = "@mimalloc",
    deps = [
        "//core:certificate_cc_proto",
        "//core:proto_io",
        "@gflags",
        "@ng-log//:ng-log",
    ],
)
