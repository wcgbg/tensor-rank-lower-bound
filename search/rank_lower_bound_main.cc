// Run the dynamic-programming search that improves the rank lower bounds in a
// certificate produced by orbit_enumerator_main.
//
// Usage:
//   bazel run --config=opt //search:rank_lower_bound_main --
//   /abs/path/orbits.pb.txt
//
// Reads a Certificate, builds an OrbitMap for the chosen problem's symmetry
// group, and runs ProcessOrbits. The result (with improved rank_lower_bound and
// rank_lower_bound_proof fields) is written to --output_path, defaulting to a
// "_updated" sibling of the input. Mirrors
// rank_search/rank_lower_bound_computer_main.cc.

#include <filesystem>
#include <limits>
#include <string>

#include <gflags/gflags.h>
#include <ng-log/logging.h>

#include "core/certificate.pb.h"
#include "core/proto_io.h"
#include "problems/chosen_problem.h"
#include "search/rank_lower_bound_computer.h"

DEFINE_string(output_path, "",
              "Output path (.pb.txt / .pb). Defaults to the input path "
              "(overwrites the input)");
DEFINE_bool(ignore_rank_lower_bound, false,
            "Clear all existing rank lower bounds before searching");
DEFINE_bool(
    ignore_non_basic_rank_lower_bound, false,
    "Clear rank lower bounds whose proof is neither flatten_matrix_proof "
    "nor forced_product_proof before searching");
DEFINE_uint64(backtracking_step_limit, std::numeric_limits<uint64_t>::max(),
              "Max backtracking search steps (all threads) per orbit; 0 to "
              "disable backtracking");
DEFINE_uint64(backtracking_max_map_size, 10'000'000,
              "Max size of the per-thread backtracking cache before stochastic "
              "halving");
DEFINE_bool(basic_method, true, "Enable the Flatten + Forced Product methods");
DEFINE_bool(degenerate_method, true, "Enable the Degenerate Reduction method");
DEFINE_int32(dim_min, 0, "Min subspace dimension to process (0 to NA)");
DEFINE_int32(dim_max, std::numeric_limits<int32_t>::max(),
             "Max subspace dimension to process (0 to NA)");
DEFINE_int32(forced_product_max_iterations_log2, 24,
             "Skip the Forced Product technique when its enumeration size "
             "num_iterations exceeds 2^this; larger enumerates more "
             "combinations (slower, potentially stronger bound)");

int main(int argc, char **argv) {
  FLAGS_alsologtostderr = true;
  FLAGS_log_dir = "/tmp";
  google::ParseCommandLineFlags(&argc, &argv, true);
  nglog::InitializeLogging(argv[0]);

  if (argc != 2) {
    LOG(FATAL) << "Usage: " << argv[0] << " <certificate_path>";
  }
  const std::string path = argv[1];

  using Problem = ChosenProblem;

  pb::Certificate certificate = ReadProtoFromFile<pb::Certificate>(path);
  CHECK_EQ(certificate.characteristic(), Problem::kP);
  CHECK_EQ(certificate.extension_degree(), Problem::kM);
  CHECK_EQ(certificate.na(), Problem::kNA);
  CHECK_EQ(certificate.nb(), Problem::kNB);
  CHECK_EQ(certificate.nc(), Problem::kNC);
  CHECK_EQ(certificate.problem_name(), Problem::Name())
      << "Certificate problem_name does not match the compiled-in problem";
  CHECK_GT(certificate.constrained_tensors_size(), 0);

  if (FLAGS_ignore_rank_lower_bound) {
    for (pb::ConstrainedTensor &rt :
         *certificate.mutable_constrained_tensors()) {
      rt.clear_rank_lower_bound();
      rt.clear_rank_lower_bound_proof();
    }
    // also remove the backtracking-proof archive
    const std::string archive_path = GetBacktrackingProofArchivePath(path);
    if (std::filesystem::exists(archive_path)) {
      std::filesystem::remove(archive_path);
    }
  }

  if (FLAGS_ignore_non_basic_rank_lower_bound) {
    for (pb::ConstrainedTensor &rt :
         *certificate.mutable_constrained_tensors()) {
      const pb::RankLowerBoundProof &proof = rt.rank_lower_bound_proof();
      if (!proof.has_flatten_matrix_proof() &&
          !proof.has_forced_product_proof()) {
        rt.set_rank_lower_bound(0);
        rt.clear_rank_lower_bound_proof();
      }
    }
  }

  std::string output_path = FLAGS_output_path;
  if (output_path.empty()) {
    output_path = path; // overwrite the input
  }

  ProcessOrbits<Problem>(
      {
          .basic_method = FLAGS_basic_method,
          .degenerate_method = FLAGS_degenerate_method,
          .backtracking_step_limit = FLAGS_backtracking_step_limit,
          .backtracking_max_map_size = FLAGS_backtracking_max_map_size,
          .dim_min = FLAGS_dim_min,
          .dim_max = FLAGS_dim_max,
          .forced_product_max_iterations_log2 =
              FLAGS_forced_product_max_iterations_log2,
      },
      output_path, &certificate);

  CHECK_GT(certificate.constrained_tensors_size(), 0);
  const pb::ConstrainedTensor &last_rt = certificate.constrained_tensors(
      certificate.constrained_tensors_size() - 1);
  CHECK(last_rt.constraints().empty());
  LOG(INFO) << "UNCONSTRAINED TENSOR RANK LOWER BOUND: "
            << last_rt.rank_lower_bound();

  LOG(INFO) << "Done. Wrote " << output_path;
  return 0;
}
