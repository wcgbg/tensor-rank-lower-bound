#include "verifier/verifier.h"

#include <string>
#include <utility>

#include "tbb/global_control.h"
#include "gtest/gtest.h"

#include "core/backtracking_proof.h"
#include "core/certificate.pb.h"
#include "core/gf_vec.h"
#include "problems/cyclic/problem.h"
#include "problems/extfield/problem.h"
#include "problems/full/problem.h"
#include "problems/matrix/problem.h"
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

int FirstBacktrackingOrbit(const pb::Certificate &certificate) {
  for (int i = 0; i < certificate.constrained_tensors_size(); ++i) {
    if (certificate.constrained_tensors(i)
            .rank_lower_bound_proof()
            .has_backtracking_proof()) {
      return i;
    }
  }
  return -1;
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

TEST(VerifierDeathTest, RejectsBacktrackingLeafPastDepthBound) {
  using Problem = full::Problem<2, 1, 3>;
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

  auto [certificate, archive_path] = BuildVerifiedCertificate<Problem>();
  const int orbit = FirstBacktrackingOrbit(certificate);
  ASSERT_GE(orbit, 0);

  BacktrackingProofArchive archive =
      BacktrackingProofArchive::Load(archive_path);
  BacktrackingProof proof = archive.Get(orbit);
  ASSERT_FALSE(proof.Empty());
  const int claimed_bound =
      certificate.constrained_tensors(orbit).rank_lower_bound();
  ASSERT_GT(claimed_bound, 0);
  ASSERT_LE(claimed_bound, 31);
  proof.dfs_constraints_size_array[0] = static_cast<uint8_t>(claimed_bound);
  archive.Set(orbit, proof);
  archive.Save(archive_path);

  EXPECT_DEATH(VerifyRankLowerBound<Problem>(certificate, archive_path),
               "exceeds the L-1 decomposition depth");
}

TEST(VerifierDeathTest, RejectsBacktrackingMaskOutsidePath) {
  using Problem = full::Problem<2, 1, 3>;
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

  auto [certificate, archive_path] = BuildVerifiedCertificate<Problem>();
  const int orbit = FirstBacktrackingOrbit(certificate);
  ASSERT_GE(orbit, 0);

  BacktrackingProofArchive archive =
      BacktrackingProofArchive::Load(archive_path);
  BacktrackingProof proof = archive.Get(orbit);
  ASSERT_FALSE(proof.Empty());
  const uint8_t depth = proof.dfs_constraints_size_array[0];
  ASSERT_GT(depth, 0);
  ASSERT_LT(depth, 32);
  proof.mask_array[0] |= uint32_t{1} << depth;
  archive.Set(orbit, proof);
  archive.Save(archive_path);

  EXPECT_DEATH(VerifyRankLowerBound<Problem>(certificate, archive_path),
               "outside the recorded path");
}

// CanonicalFromWitness must reject witnesses that are not genuine group
// elements as encoded in the certificate, before applying them.
TEST(VerifierDeathTest, RejectsMalformedWitness) {
  using Matrix = matrix::Problem<2, 1, 2, 2, 3>; // non-cubic: no transpose
  const Matrix::SymmetryGroup matrix_group;
  const Constraints<2, 1, Matrix::kNA> q = {DecodeGFVec<2, 1, Matrix::kNA>(1)};
  const uint32_t query_id =
      static_cast<uint32_t>(matrix_group.query.Identity());
  const uint32_t store_id = matrix_group.store.Identity();
  EXPECT_EQ(CanonicalFromWitness<Matrix>(matrix_group, q, query_id, store_id),
            q);
  // High bits that a uint16 store element would silently drop.
  EXPECT_DEATH(CanonicalFromWitness<Matrix>(matrix_group, q, query_id,
                                            store_id | 0x10000u),
               "does not fit its encoding");
  // Singular R / L.
  EXPECT_DEATH(CanonicalFromWitness<Matrix>(matrix_group, q, query_id, 0),
               "store element is not in the symmetry group");
  EXPECT_DEATH(CanonicalFromWitness<Matrix>(matrix_group, q, 0, store_id),
               "query element is not in the symmetry group");
  // Transpose flag on a non-cubic format.
  EXPECT_DEATH(CanonicalFromWitness<Matrix>(matrix_group, q,
                                            query_id | (1u << 16), store_id),
               "query element is not in the symmetry group");

  using Full = full::Problem<2, 1, 3>;
  const Full::SymmetryGroup full_group;
  const Constraints<2, 1, Full::kNA> fq = {DecodeGFVec<2, 1, Full::kNA>(1)};
  EXPECT_DEATH(CanonicalFromWitness<Full>(
                   full_group, fq, 0,
                   static_cast<uint32_t>(full_group.store.Size())),
               "store element is not in the symmetry group");
  EXPECT_DEATH(CanonicalFromWitness<Full>(full_group, fq, 0, 0x80000000u),
               "store element is not in the symmetry group");
}

// Forced product on a tensor whose search space (q^(r2p * r1)) overflows
// uint64 must not report a bound. full N=12 restricted to A-support {0, 6} has,
// in the C-mode, 12 independent rank-1 slices and 6 others: 2^72 combinations.
// Its rank is at most 2 * 12 (two nonzero A-slices).
TEST(VerifierDeathTest, ForcedProductRejectsOverflowingSearch) {
  using Problem = full::Problem<2, 1, 12>;
  constexpr int NA = Problem::kNA;
  constexpr int NB = Problem::kNB;
  constexpr int NC = Problem::kNC;
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

  Constraints<2, 1, NA> constraints;
  for (int i = 0; i < NA; ++i) {
    if (i != 0 && i != 6) {
      constraints.push_back(DecodeGFVec<2, 1, NA>(uint64_t{1} << i));
    }
  }
  GaussJordanRREF<2, 1, NA>(&constraints);
  const Tensor<2, 1, NA, NB, NC> tensor =
      ApplyConstraintsToTensor<2, 1, NA, NB, NC>(constraints,
                                                 Problem::MakeTensor());
  const Tensor<2, 1, NC, NA, NB> c_mode = CyclicTranspose<2, 1, NB, NC, NA>(
      CyclicTranspose<2, 1, NA, NB, NC>(tensor));
  EXPECT_LE((RankLowerBoundForcedProductA<2, 1, NC, NA, NB>(c_mode)), 2 * NB);

  pb::ForcedProductProof proof;
  proof.set_projection_type(2);
  EXPECT_DEATH((VerifyForcedProductProof<2, 1, NA, NB, NC>(tensor, 1000, proof)),
               "");
}

// Certificate rows must be elements of F_q^NA with no zero rows, and a
// degenerate proof's extra constraint must be in range; otherwise the checks
// would index tensors and field tables out of range.
TEST(VerifierDeathTest, RejectsMalformedConstraints) {
  using Problem = full::Problem<2, 1, 3>; // NA = 3, rows are uint8
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

  const auto [certificate, archive_path] = BuildVerifiedCertificate<Problem>();
  int with_constraints = -1;
  for (int i = 0; i < certificate.constrained_tensors_size(); ++i) {
    const size_t rows = certificate.constrained_tensors(i).constraints().size();
    if (rows > 0 && rows < Problem::kNA) { // room for one more row
      with_constraints = i;
      break;
    }
  }
  ASSERT_GE(with_constraints, 0);

  {
    pb::Certificate bad = certificate;
    std::string *bytes =
        bad.mutable_constrained_tensors(with_constraints)->mutable_constraints();
    (*bytes)[0] = static_cast<char>((*bytes)[0] | 0x80); // bit 7 >= NA
    EXPECT_DEATH(VerifyRankLowerBound<Problem>(bad, archive_path),
                 "outside F_q\\^NA");
  }
  {
    pb::Certificate bad = certificate;
    std::string *bytes =
        bad.mutable_constrained_tensors(with_constraints)->mutable_constraints();
    bytes->insert(bytes->begin(), '\0'); // a leading zero row is still RREF
    EXPECT_DEATH(VerifyRankLowerBound<Problem>(bad, archive_path),
                 "zero row");
  }
  {
    pb::Certificate bad = certificate;
    pb::DegenerateProof *proof =
        bad.mutable_constrained_tensors(with_constraints)
            ->mutable_rank_lower_bound_proof()
            ->mutable_degenerate_proof();
    proof->set_extra_constraint(0x81); // bit 7 >= NA
    EXPECT_DEATH(VerifyRankLowerBound<Problem>(bad, archive_path),
                 "extra_constraint is not an element");
  }
}

// A certificate must name the problem it was built for; an empty or different
// name is rejected even when the dimensions match.
TEST(VerifierDeathTest, RequiresMatchingProblemName) {
  using Problem = cyclic::Problem<3, 1, 3>;
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

  auto [certificate, archive_path] = BuildVerifiedCertificate<Problem>();
  certificate.clear_problem_name();
  EXPECT_DEATH(VerifyRankLowerBound<Problem>(certificate, archive_path),
               "problem_name does not match");
}

} // namespace
