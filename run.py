#!/usr/bin/env python3
"""Drive the four-stage prover pipeline for a chosen problem.

The first positional argument selects which stages run:
  search-all             orbit_enumerator_main, rank_upper_bound_main,
                         rank_lower_bound_main for every problem listed in this file
                         (please uncomment the problems you want to run).
  search-and-verify-all  search-all then verify, for every problem listed in this file.
  verify <cert.pb.txt>   verifier_main on the single certificate <cert.pb.txt>; the
                         problem and its parameters are derived from the filename
"""

import argparse
import logging
import os
import subprocess
import time

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
assert os.path.exists(os.path.join(REPO_ROOT, "MODULE.bazel")), REPO_ROOT
RUN_LOG = os.path.join(REPO_ROOT, "log", "run.log")


def chosen_problem_copt(defines):
    # Select the compiled-in ChosenProblem at build time via -D defines (see
    # chosen_problem.h). chosen_problem.h is included only by the *_main.cc
    # binaries, so scope the defines to those translation units to avoid
    # invalidating the rest of the build cache. The defines are comma-free
    # individually; the commas here separate distinct copts (Bazel splits
    # --per_file_copt option lists on commas).
    opts = ",".join(f"-D{d}" for d in defines)
    return f"--per_file_copt=.*_main\\.cc@{opts}"


dry_run = False
command = None

log = logging.getLogger(__name__)


def configure_logging():
    if log.handlers:
        return
    os.makedirs(os.path.dirname(RUN_LOG), exist_ok=True)
    log.setLevel(logging.INFO)
    handler = logging.FileHandler(RUN_LOG)
    handler.setFormatter(
        logging.Formatter("%(asctime)s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    )
    log.addHandler(handler)


def sh(cmd):
    cmd_str = " ".join(cmd)
    print("+ " + cmd_str, flush=True)
    log.info("RUN     %s", cmd_str)
    if dry_run:
        return 0.0
    t0 = time.perf_counter()
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)
    return time.perf_counter() - t0


def search_pipeline(problem, defines, q, lb_args_list, ub_args, cert=None):
    # Builds the binaries needed for the selected command (with the chosen
    # problem injected via -DCP_* defines) and runs the corresponding stages.
    # "search-all" runs enumerate -> (upper) -> lower; "verify" runs verify;
    # "search-and-verify-all" runs both. When cert is None the certificate path
    # is computed as certs/<family>/cert_<problem>.pb.txt, with <family> the
    # leading token of the problem name; the "verify" command passes an explicit
    # cert path instead.
    family = problem.split("_", 1)[0]
    if cert is None:
        cert = f"certs/{family}/cert_{problem}.pb.txt"

    do_search = command in ("search-all", "search-and-verify-all")
    do_verify = command in ("verify", "search-and-verify-all")

    build_targets = []
    if do_search:
        build_targets += [
            "//search:orbit_enumerator_main",
            "//search:rank_lower_bound_main",
        ]
        if q == 2:
            build_targets.append("//upper_bound:rank_upper_bound_main")
    if do_verify:
        build_targets.append("//verifier:verifier_main")
    sh(["bazel", "build", "--config=opt", chosen_problem_copt(defines)] + build_targets)

    search_time = 0.0
    if do_search:
        search_time += sh(
            ["bazel-bin/search/orbit_enumerator_main", "--output_path", cert]
        )

        # The flip-graph upper-bound search implementation is F_2-only.
        if q == 2:
            if ub_args is not None:
                cmd = ["bazel-bin/upper_bound/rank_upper_bound_main", cert]
                cmd.extend([f"--{k}={v}" for k, v in ub_args.items()])
                sh(cmd)
        else:
            assert ub_args is None

        for i, lb_args in enumerate(lb_args_list):
            lower_cmd = ["bazel-bin/search/rank_lower_bound_main", cert]
            assert lb_args is not None
            lower_cmd.extend([f"--{k}={v}" for k, v in lb_args.items()])
            if i == 0:
                lower_cmd.append("--ignore_rank_lower_bound")
            search_time += sh(lower_cmd)

    verification_time = 0.0
    if do_verify:
        verification_time = sh(["bazel-bin/verifier/verifier_main", cert])

    log.info(
        "RESULT  %s  search_time=%.3f s  verification_time=%.3f s",
        problem,
        search_time,
        verification_time,
    )


def search_polynomial(family, p, m, n, lb_args, ub_args):
    q = p**m
    problem = f"{family}_q{q:02d}_n{n:02d}"
    defines = [f"CP_{family.upper()}", f"CP_P={p}", f"CP_M={m}", f"CP_N={n}"]

    print()
    print("========", problem, "========")

    search_pipeline(problem, defines, q, [lb_args], ub_args)


def search_matrix(p, m, n0, n1, n2, lb_args_pass0, lb_args_pass1, ub_args):
    q = p**m
    problem = f"matrix_q{q:02d}_n{n0}{n1}{n2}"
    defines = [
        "CP_MATRIX",
        f"CP_P={p}",
        f"CP_M={m}",
        f"CP_N0={n0}",
        f"CP_N1={n1}",
        f"CP_N2={n2}",
    ]

    print()
    print("========", problem, "========")

    assert lb_args_pass0 is not None
    lb_args_list = [lb_args_pass0]
    if lb_args_pass1 is not None:
        lb_args_list.append(lb_args_pass1)
    search_pipeline(problem, defines, q, lb_args_list, ub_args)


def factor_prime_power(q):
    # Return (p, m) with q == p**m. The smallest divisor >= 2 of q is always
    # prime, and a prime power has a unique prime base, so (p, m) is uniquely
    # determined (see summary/prev_results_polynomial/gen_prev_results.py:prime_powers).
    p = next(d for d in range(2, q + 1) if q % d == 0)
    t, m = q, 0
    while t % p == 0:
        t //= p
        m += 1
    assert t == 1, f"{q} is not a prime power"
    return p, m


def parse_cert(cert_path):
    # Derive (problem, defines, q) from a certificate filename such as
    # certs/matrix/cert_matrix_q02_n333.pb.txt -> problem "matrix_q02_n333".
    # The build-time -DCP_* problem-selection defines mirror those built by
    # search_polynomial()/search_matrix().
    name = os.path.basename(cert_path)
    assert name.startswith("cert_") and name.endswith(".pb.txt"), cert_path
    problem = name[len("cert_") : -len(".pb.txt")]

    tokens = problem.split("_")
    family = tokens[0]
    q_token, n_token = tokens[1], tokens[2]
    assert q_token.startswith("q") and n_token.startswith("n"), problem
    q = int(q_token[1:])
    p, m = factor_prime_power(q)
    defines = [f"CP_P={p}", f"CP_M={m}"]

    if family == "matrix":
        digits = n_token[1:]
        # The naming scheme f"matrix_q{q:02d}_n{n0}{n1}{n2}" concatenates the
        # three single-digit matrix sizes.
        assert len(digits) == 3, problem
        n0, n1, n2 = (int(d) for d in digits)
        defines = ["CP_MATRIX", *defines, f"CP_N0={n0}", f"CP_N1={n1}", f"CP_N2={n2}"]
    else:
        n = int(n_token[1:])
        defines = [f"CP_{family.upper()}", *defines, f"CP_N={n}"]

    return problem, defines, q


def verify_cert(cert_path):
    if not dry_run:
        assert os.path.exists(cert_path), cert_path
    problem, defines, q = parse_cert(cert_path)

    print()
    print("========", problem, "========")

    search_pipeline(problem, defines, q, [], None, cert=cert_path)


def search_full_all():
    pass
    # # q=2
    # search_polynomial("full", 2, 1, 1, {}, {})
    # search_polynomial("full", 2, 1, 2, {}, {})
    # search_polynomial("full", 2, 1, 3, {}, {})
    # search_polynomial("full", 2, 1, 4, {}, {})
    # lb_args = {"backtracking_step_limit": 1_000_000}
    # search_polynomial("full", 2, 1, 5, lb_args, {})
    # lb_args = {
    #     "backtracking_step_limit": 1_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # search_polynomial("full", 2, 1, 6, lb_args, {})
    # lb_args = {
    #     "backtracking_step_limit": 1_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 1_000_000,
    #     "max_steps_at_a_rank": 1_000,
    # }
    # search_polynomial("full", 2, 1, 7, lb_args, ub_args)
    # lb_args = {
    #     "backtracking_step_limit": 1_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 10_000,
    # }
    # search_polynomial("full", 2, 1, 8, lb_args, ub_args)
    # # q=3
    # search_polynomial("full", 3, 1, 1, {}, None)
    # search_polynomial("full", 3, 1, 2, {}, None)
    # search_polynomial("full", 3, 1, 3, {}, None)
    # search_polynomial("full", 3, 1, 4, {}, None)
    # lb_args = {
    #     "backtracking_step_limit": 1_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # search_polynomial("full", 3, 1, 5, lb_args, None)
    # lb_args = {
    #     "backtracking_step_limit": 1_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # search_polynomial("full", 3, 1, 6, lb_args, None)
    # # q=4
    # search_polynomial("full", 2, 2, 1, {}, None)
    # search_polynomial("full", 2, 2, 2, {}, None)
    # search_polynomial("full", 2, 2, 3, {}, None)
    # search_polynomial("full", 2, 2, 4, {}, None)
    # # q=5
    # search_polynomial("full", 5, 1, 1, {}, None)
    # search_polynomial("full", 5, 1, 2, {}, None)
    # search_polynomial("full", 5, 1, 3, {}, None)
    # # q=7
    # search_polynomial("full", 7, 1, 1, {}, None)
    # search_polynomial("full", 7, 1, 2, {}, None)
    # search_polynomial("full", 7, 1, 3, {}, None)
    # # q=8
    # search_polynomial("full", 2, 3, 1, {}, None)
    # search_polynomial("full", 2, 3, 2, {}, None)
    # search_polynomial("full", 2, 3, 3, {}, None)
    # search_polynomial("full", 2, 3, 4, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("full", 2, 3, 5, lb_args, None)
    # # q=9
    # search_polynomial("full", 3, 2, 1, {}, None)
    # search_polynomial("full", 3, 2, 2, {}, None)
    # search_polynomial("full", 3, 2, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("full", 3, 2, 4, lb_args, None)
    # # q=11
    # search_polynomial("full", 11, 1, 1, {}, None)
    # search_polynomial("full", 11, 1, 2, {}, None)
    # search_polynomial("full", 11, 1, 3, {}, None)
    # # q=13
    # search_polynomial("full", 13, 1, 1, {}, None)
    # search_polynomial("full", 13, 1, 2, {}, None)
    # search_polynomial("full", 13, 1, 3, {}, None)
    # # q=16
    # search_polynomial("full", 2, 4, 1, {}, None)
    # search_polynomial("full", 2, 4, 2, {}, None)
    # search_polynomial("full", 2, 4, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("full", 2, 4, 4, lb_args, None)


def search_cyclic_all():
    pass
    # # q=2
    # search_polynomial("cyclic", 2, 1, 1, {}, {})
    # search_polynomial("cyclic", 2, 1, 2, {}, {})
    # search_polynomial("cyclic", 2, 1, 3, {}, {})
    # search_polynomial("cyclic", 2, 1, 4, {}, {})
    # search_polynomial("cyclic", 2, 1, 5, {}, {})
    # search_polynomial("cyclic", 2, 1, 6, {}, {})
    # search_polynomial("cyclic", 2, 1, 7, {}, {})
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 10_000,
    #     "num_paths": 192,
    # }
    # search_polynomial("cyclic", 2, 1, 8, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 100_000_000,
    #     "max_steps_at_a_rank": 10_000,
    #     "num_paths": 192,
    # }
    # search_polynomial("cyclic", 2, 1, 9, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 100_000_000,
    #     "max_steps_at_a_rank": 10_000,
    #     "num_paths": 192,
    # }
    # search_polynomial("cyclic", 2, 1, 10, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 1, 11, lb_args, None)
    # # q=3
    # search_polynomial("cyclic", 3, 1, 1, {}, None)
    # search_polynomial("cyclic", 3, 1, 2, {}, None)
    # search_polynomial("cyclic", 3, 1, 3, {}, None)
    # search_polynomial("cyclic", 3, 1, 4, {}, None)
    # search_polynomial("cyclic", 3, 1, 5, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 3, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 3, 1, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 3, 1, 8, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 3, 1, 9, lb_args, None)
    # # q=4
    # search_polynomial("cyclic", 2, 2, 1, {}, None)
    # search_polynomial("cyclic", 2, 2, 2, {}, None)
    # search_polynomial("cyclic", 2, 2, 3, {}, None)
    # search_polynomial("cyclic", 2, 2, 4, {}, None)
    # search_polynomial("cyclic", 2, 2, 5, {}, None)
    # search_polynomial("cyclic", 2, 2, 6, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 2, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 2, 8, lb_args, None)
    # # q=5
    # search_polynomial("cyclic", 5, 1, 1, {}, None)
    # search_polynomial("cyclic", 5, 1, 2, {}, None)
    # search_polynomial("cyclic", 5, 1, 3, {}, None)
    # search_polynomial("cyclic", 5, 1, 4, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 5, 1, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 5, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 5, 1, 7, lb_args, None)
    # # q=7
    # search_polynomial("cyclic", 7, 1, 1, {}, None)
    # search_polynomial("cyclic", 7, 1, 2, {}, None)
    # search_polynomial("cyclic", 7, 1, 3, {}, None)
    # search_polynomial("cyclic", 7, 1, 4, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 7, 1, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 7, 1, 6, lb_args, None)
    # # q=8
    # search_polynomial("cyclic", 2, 3, 1, {}, None)
    # search_polynomial("cyclic", 2, 3, 2, {}, None)
    # search_polynomial("cyclic", 2, 3, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 3, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 3, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 3, 6, lb_args, None)
    # # q=9
    # search_polynomial("cyclic", 3, 2, 1, {}, None)
    # search_polynomial("cyclic", 3, 2, 2, {}, None)
    # search_polynomial("cyclic", 3, 2, 3, {}, None)
    # search_polynomial("cyclic", 3, 2, 4, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 3, 2, 5, lb_args, None)
    # # q=11
    # search_polynomial("cyclic", 11, 1, 1, {}, None)
    # search_polynomial("cyclic", 11, 1, 2, {}, None)
    # search_polynomial("cyclic", 11, 1, 3, {}, None)
    # search_polynomial("cyclic", 11, 1, 4, {}, None)
    # search_polynomial("cyclic", 11, 1, 5, {}, None)
    # # q=13
    # search_polynomial("cyclic", 13, 1, 1, {}, None)
    # search_polynomial("cyclic", 13, 1, 2, {}, None)
    # search_polynomial("cyclic", 13, 1, 3, {}, None)
    # search_polynomial("cyclic", 13, 1, 4, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 13, 1, 5, lb_args, None)
    # # q=16
    # search_polynomial("cyclic", 2, 4, 1, {}, None)
    # search_polynomial("cyclic", 2, 4, 2, {}, None)
    # search_polynomial("cyclic", 2, 4, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 4, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 4, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("cyclic", 2, 4, 6, lb_args, None)


def search_extfield_all():
    pass
    # # q=2
    # search_polynomial("extfield", 2, 1, 1, {}, {})
    # search_polynomial("extfield", 2, 1, 2, {}, {})
    # search_polynomial("extfield", 2, 1, 3, {}, {})
    # search_polynomial("extfield", 2, 1, 4, {}, {})
    # lb_args = {"backtracking_step_limit": 10_000_000}
    # search_polynomial("extfield", 2, 1, 5, lb_args, {})
    # lb_args = {"backtracking_step_limit": 10_000_000}
    # search_polynomial("extfield", 2, 1, 6, lb_args, {})
    # lb_args = {"backtracking_step_limit": 10_000_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 1_000,
    # }
    # search_polynomial("extfield", 2, 1, 7, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 1, 8, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 1, 9, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 1, 10, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 1, 11, lb_args, None)
    # # q=3
    # search_polynomial("extfield", 3, 1, 1, {}, None)
    # search_polynomial("extfield", 3, 1, 2, {}, None)
    # search_polynomial("extfield", 3, 1, 3, {}, None)
    # search_polynomial("extfield", 3, 1, 4, {}, None)
    # search_polynomial("extfield", 3, 1, 5, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 3, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 3, 1, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 3, 1, 8, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 3, 1, 9, lb_args, None)
    # # q=4
    # search_polynomial("extfield", 2, 2, 1, {}, None)
    # search_polynomial("extfield", 2, 2, 2, {}, None)
    # search_polynomial("extfield", 2, 2, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 2, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 2, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 2, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 2, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 2, 8, lb_args, None)
    # # q=5
    # search_polynomial("extfield", 5, 1, 1, {}, None)
    # search_polynomial("extfield", 5, 1, 2, {}, None)
    # search_polynomial("extfield", 5, 1, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 5, 1, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 5, 1, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 5, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 5, 1, 7, lb_args, None)
    # # q=7
    # search_polynomial("extfield", 7, 1, 1, {}, None)
    # search_polynomial("extfield", 7, 1, 2, {}, None)
    # search_polynomial("extfield", 7, 1, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 7, 1, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 7, 1, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 7, 1, 6, lb_args, None)
    # # q=8
    # search_polynomial("extfield", 2, 3, 1, {}, None)
    # search_polynomial("extfield", 2, 3, 2, {}, None)
    # search_polynomial("extfield", 2, 3, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 3, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 3, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 3, 6, lb_args, None)
    # # q=9
    # search_polynomial("extfield", 3, 2, 1, {}, None)
    # search_polynomial("extfield", 3, 2, 2, {}, None)
    # search_polynomial("extfield", 3, 2, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 3, 2, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 3, 2, 5, lb_args, None)
    # # q=11
    # search_polynomial("extfield", 11, 1, 1, {}, None)
    # search_polynomial("extfield", 11, 1, 2, {}, None)
    # search_polynomial("extfield", 11, 1, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 11, 1, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 11, 1, 5, lb_args, None)
    # # q=13
    # search_polynomial("extfield", 13, 1, 1, {}, None)
    # search_polynomial("extfield", 13, 1, 2, {}, None)
    # search_polynomial("extfield", 13, 1, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 13, 1, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 13, 1, 5, lb_args, None)
    # # q=16
    # search_polynomial("extfield", 2, 4, 1, {}, None)
    # search_polynomial("extfield", 2, 4, 2, {}, None)
    # search_polynomial("extfield", 2, 4, 3, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 4, 4, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("extfield", 2, 4, 5, lb_args, None)


def search_matrix_all():
    pass
    # # q=2
    # search_matrix(2, 1, 2, 2, 2, {}, None, {})
    # search_matrix(2, 1, 2, 2, 3, {}, None, {})
    # search_matrix(2, 1, 2, 2, 4, {}, None, {})
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_matrix(2, 1, 2, 3, 3, lb_args, None, {})
    # lb_args = {
    #     "backtracking_step_limit": 10_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # search_matrix(2, 1, 3, 2, 4, lb_args, None, {})
    # lb_args = {
    #     "backtracking_step_limit": 10_000_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 1_000_000,
    #     "max_steps_at_a_rank": 1_000,
    #     "num_paths": 8,
    # }
    # search_matrix(2, 1, 3, 3, 3, lb_args, None, ub_args)
    # lb_args_pass0 = {
    #     "backtracking_step_limit": 100_000_000,
    #     "forced_product_max_iterations_log2": 32,
    #     "dim_min": 4,
    # }
    # lb_args_pass1 = {
    #     "backtracking_step_limit": 100_000_000,
    #     "forced_product_max_iterations_log2": 32,
    #     "dim_max": 3,
    #     "backtracking_max_map_size": 3_000_000,  # reduce memory
    # }
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 10_000,
    #     "num_paths": 16,
    # }
    # search_matrix(2, 1, 3, 3, 4, lb_args_pass0, lb_args_pass1, ub_args)
    # lb_args = {
    #     "backtracking_step_limit": 100_000,
    #     "forced_product_max_iterations_log2": 32,
    # }
    # search_matrix(2, 1, 3, 4, 4, lb_args, None, None)


def search_truncated_all():
    # truncated = short product, F_q[x]/(x^n)
    pass
    # # q=2
    # search_polynomial("truncated", 2, 1, 1, {}, {})
    # search_polynomial("truncated", 2, 1, 2, {}, {})
    # search_polynomial("truncated", 2, 1, 3, {}, {})
    # search_polynomial("truncated", 2, 1, 4, {}, {})
    # search_polynomial("truncated", 2, 1, 5, {}, {})
    # search_polynomial("truncated", 2, 1, 6, {}, {})
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 2, 1, 7, lb_args, {})
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 1_000_000,
    #     "max_steps_at_a_rank": 1_000,
    # }
    # search_polynomial("truncated", 2, 1, 8, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 1_000_000,
    #     "max_steps_at_a_rank": 1_000,
    # }
    # search_polynomial("truncated", 2, 1, 9, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 10_000,
    # }
    # search_polynomial("truncated", 2, 1, 10, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 10_000,
    #     "num_paths": 192,
    # }
    # search_polynomial("truncated", 2, 1, 11, lb_args, ub_args)
    # lb_args = {"backtracking_step_limit": 100_000}
    # ub_args = {
    #     "last_only": "true",
    #     "path_limit": 10_000_000,
    #     "max_steps_at_a_rank": 10_000,
    #     "num_paths": 192,
    # }
    # search_polynomial("truncated", 2, 1, 12, lb_args, ub_args)
    # # q=3
    # search_polynomial("truncated", 3, 1, 1, {}, None)
    # search_polynomial("truncated", 3, 1, 2, {}, None)
    # search_polynomial("truncated", 3, 1, 3, {}, None)
    # search_polynomial("truncated", 3, 1, 4, {}, None)
    # search_polynomial("truncated", 3, 1, 5, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 3, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 3, 1, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 3, 1, 8, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 3, 1, 9, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 3, 1, 10, lb_args, None)
    # # q=4
    # search_polynomial("truncated", 2, 2, 1, {}, None)
    # search_polynomial("truncated", 2, 2, 2, {}, None)
    # search_polynomial("truncated", 2, 2, 3, {}, None)
    # search_polynomial("truncated", 2, 2, 4, {}, None)
    # search_polynomial("truncated", 2, 2, 5, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 2, 2, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 2, 2, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 2, 2, 8, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 2, 2, 9, lb_args, None)
    # # q=5
    # search_polynomial("truncated", 5, 1, 1, {}, None)
    # search_polynomial("truncated", 5, 1, 2, {}, None)
    # search_polynomial("truncated", 5, 1, 3, {}, None)
    # search_polynomial("truncated", 5, 1, 4, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 5, 1, 5, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 5, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 5, 1, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("truncated", 5, 1, 8, lb_args, None)


def search_negacyclic_all():
    # negacyclic, F_q[x]/(x^n+1); over char 2 it equals cyclic
    pass
    # # q=3
    # search_polynomial("negacyclic", 3, 1, 1, {}, None)
    # search_polynomial("negacyclic", 3, 1, 2, {}, None)
    # search_polynomial("negacyclic", 3, 1, 3, {}, None)
    # search_polynomial("negacyclic", 3, 1, 4, {}, None)
    # search_polynomial("negacyclic", 3, 1, 5, {}, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("negacyclic", 3, 1, 6, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("negacyclic", 3, 1, 7, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("negacyclic", 3, 1, 8, lb_args, None)
    # lb_args = {"backtracking_step_limit": 100_000}
    # search_polynomial("negacyclic", 3, 1, 9, lb_args, None)


def main():
    global dry_run, command

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "command",
        choices=["search-all", "verify", "search-and-verify-all"],
    )
    parser.add_argument(
        "cert",
        nargs="?",
        help="certificate to verify (required for the 'verify' command)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print commands without running them",
    )
    args = parser.parse_args()
    dry_run = args.dry_run
    command = args.command
    configure_logging()

    if command == "verify":
        if not args.cert:
            parser.error("the 'verify' command requires a certificate path argument")
        verify_cert(args.cert)
        return

    if args.cert:
        parser.error(f"the '{command}' command takes no positional cert argument")

    search_full_all()
    search_cyclic_all()
    search_extfield_all()
    search_matrix_all()
    search_truncated_all()
    search_negacyclic_all()


if __name__ == "__main__":
    main()
