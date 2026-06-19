#include "verifier/verifier.h"

#include <string>
#include <utility>

#include <gtest/gtest.h>
#include <tbb/global_control.h>

#include "core/backtracking_proof.h"
#include "core/certificate.pb.h"
#include "problems/cyclic/problem.h"
#include "problems/extfield/problem.h"
#include "problems/full/problem.h"
#include "search/orbit_enumerator_slow.h"
#include "search/rank_lower_bound_computer.h"

namespace {

// Enumerate orbits and run the full DP search, writing the certificate (and the
// backtracking-proof archive) under a temporary directory.
template <class Problem>
std::pair<pb::Certificate, std::string> BuildVerifiedCertificate() {
  const std::string temp_dir = testing::TempDir();
  const std::string output_path = temp_dir + "/" + Problem::Name() + ".pb.txt";

  const typename Problem::SymmetryGroup group;
  const OrbitEnumeratorSlow<Problem> enumerator(&group);
  pb::Certificate certificate = enumerator.Search();

  ProcessOrbits<Problem>(ProcessOptions{}, output_path, &certificate);
  return {certificate, GetBacktrackingProofArchivePath(output_path)};
}

template <class Problem> int RunPipeline() {
  auto [certificate, archive_path] = BuildVerifiedCertificate<Problem>();
  if (certificate.constrained_tensors_size() == 0)
    return -1;
  EXPECT_EQ(certificate.characteristic(), Problem::kP);
  return VerifyRankLowerBound<Problem>(certificate, archive_path);
}

TEST(VerifierTest, VerifiesFullPolyMulN3) {
  using Problem = full::Problem<2, 1, 3>;
  // n=3 full polynomial multiplication has bilinear complexity 6 over F_2.
  EXPECT_EQ(RunPipeline<Problem>(), 6);
}

TEST(VerifierTest, CyclicN3OverF3) {
  using Problem = cyclic::Problem<3, 1, 3>;
  EXPECT_GT(RunPipeline<Problem>(), 0);
}

TEST(VerifierTest, FullN2OverF3) {
  using Problem = full::Problem<3, 1, 2>;
  // Multiplying two degree-1 polynomials over F_3: rank is 3 (Karatsuba).
  EXPECT_EQ(RunPipeline<Problem>(), 3);
}

TEST(VerifierTest, ExtfieldN2OverF3) {
  using Problem = extfield::Problem<3, 1, 2>;
  // F_9 multiplication over F_3: bilinear rank ≥ 3 (Karatsuba-style is tight
  // at 3 for the 2x2x2 tensor).
  EXPECT_EQ(RunPipeline<Problem>(), 3);
}

// M ≥ 2: base field 𝔽_q is itself an extension. These exercise the full
// pipeline (orbit enumeration → all four lower-bound techniques → verification
// with the projective store inverse and the semilinear query dual) over the
// 𝔽_q-aware symmetry group.
TEST(VerifierTest, ExtfieldN2OverF4) {
  // 𝔽_16 = 𝔽_{4^2} over 𝔽_4. The degree-2 extension multiplication tensor has
  // bilinear rank 3 (the 2n−1 bound, tight for n=2) over any base field.
  using Problem = extfield::Problem<2, 2, 2>;
  EXPECT_EQ(RunPipeline<Problem>(), 3);
}

TEST(VerifierTest, ExtfieldN2OverF9) {
  // 𝔽_81 = 𝔽_{9^2} over 𝔽_9, again rank 3.
  using Problem = extfield::Problem<3, 2, 2>;
  EXPECT_EQ(RunPipeline<Problem>(), 3);
}

TEST(VerifierDeathTest, DetectsInflatedBound) {
  using Problem = full::Problem<2, 1, 3>;
  // EXPECT_DEATH only works with a single thread.
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

  auto [certificate, archive_path] = BuildVerifiedCertificate<Problem>();
  ASSERT_GT(certificate.constrained_tensors_size(), 0);
  pb::ConstrainedTensor *last = certificate.mutable_constrained_tensors(
      certificate.constrained_tensors_size() - 1);
  last->set_rank_lower_bound(last->rank_lower_bound() + 1);
  EXPECT_DEATH(VerifyRankLowerBound<Problem>(certificate, archive_path), "");
}

} // namespace
