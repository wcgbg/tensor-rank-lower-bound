// Summarize every cached certificate under a directory tree.
//
// Usage:
//   bazel run --config=opt //summary/our_results:cert_summary_main
//   bazel run --config=opt //summary/our_results:cert_summary_main -- --certs_dir=certs
//
// Scans `--certs_dir` recursively for *.pb.txt files, parses each as a
// pb::Certificate, and reports the headline fields (problem_name,
// characteristic, extension_degree, na, nb, nc) together with the
// rank_lower_bound / rank_upper_bound of the LAST ConstrainedTensor (the
// unconstrained tensor; see CLAUDE.md). The summary is written as
// cert_summary.csv and cert_summary.md inside `--out_dir`.

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <gflags/gflags.h>
#include <ng-log/logging.h>

#include "core/certificate.pb.h"
#include "core/proto_io.h"

DEFINE_string(certs_dir, "certs",
              "Directory scanned recursively for *.pb.txt certificates");
DEFINE_string(
    out_dir, "summary/our_results",
    "Directory where cert_summary.csv and cert_summary.md are written");

namespace {

namespace fs = std::filesystem;

// Resolve a possibly-relative path against the directory `bazel run` was
// invoked from (BUILD_WORKING_DIRECTORY), so the defaults work from the repo
// root. Absolute paths and runs without that env var are left unchanged.
fs::path ResolvePath(const std::string &path) {
  fs::path p(path);
  if (p.is_absolute()) {
    return p;
  }
  if (const char *cwd = std::getenv("BUILD_WORKING_DIRECTORY")) {
    return fs::path(cwd) / p;
  }
  return p;
}

struct Row {
  std::string file; // path relative to the scanned root
  std::string problem_name;
  int p = 0;
  int m = 0;
  long long q = 0; // p^m
  int na = 0;
  int nb = 0;
  int nc = 0;
  std::string lb; // "-" when unset
  std::string ub; // "-" when unset
};

long long IntPow(int base, int exp) {
  long long r = 1;
  for (int i = 0; i < exp; ++i) {
    r *= base;
  }
  return r;
}

// Column headers, shared by both outputs.
const std::vector<std::string> &Headers() {
  static const std::vector<std::string> headers = {
      "file", "problem", "p", "m", "q", "na", "nb", "nc", "lb", "ub"};
  return headers;
}

// One row of cell values, in the same column order as Headers().
std::vector<std::string> Cells(const Row &r) {
  return {r.file,
          r.problem_name,
          std::to_string(r.p),
          std::to_string(r.m),
          std::to_string(r.q),
          std::to_string(r.na),
          std::to_string(r.nb),
          std::to_string(r.nc),
          r.lb,
          r.ub};
}

void WriteCsv(const fs::path &path, const std::vector<Row> &rows) {
  std::ofstream csv(path);
  CHECK(csv) << "Failed to open CSV for writing: " << path.string();
  auto write_line = [&](const std::vector<std::string> &c) {
    for (std::size_t i = 0; i < c.size(); ++i) {
      csv << (i ? "," : "") << c[i];
    }
    csv << '\n';
  };
  write_line(Headers());
  for (const Row &r : rows) {
    write_line(Cells(r));
  }
}

// Writes a Markdown table, with columns padded to their max width so the raw
// source is readable too.
void WriteMarkdown(const fs::path &path, const std::vector<Row> &rows) {
  const std::vector<std::string> &headers = Headers();
  std::vector<std::size_t> width(headers.size());
  for (std::size_t i = 0; i < headers.size(); ++i) {
    width[i] = headers[i].size();
  }
  for (const Row &r : rows) {
    const std::vector<std::string> c = Cells(r);
    for (std::size_t i = 0; i < c.size(); ++i) {
      width[i] = std::max(width[i], c[i].size());
    }
  }

  std::ofstream md(path);
  CHECK(md) << "Failed to open Markdown for writing: " << path.string();
  auto pad = [](const std::string &s, std::size_t w) {
    return s + std::string(w - s.size(), ' ');
  };
  auto write_row = [&](const std::vector<std::string> &c) {
    md << "|";
    for (std::size_t i = 0; i < c.size(); ++i) {
      md << ' ' << pad(c[i], width[i]) << " |";
    }
    md << '\n';
  };
  write_row(headers);
  md << "|";
  for (std::size_t i = 0; i < headers.size(); ++i) {
    md << ' ' << std::string(width[i], '-') << " |";
  }
  md << '\n';
  for (const Row &r : rows) {
    write_row(Cells(r));
  }
}

} // namespace

int main(int argc, char **argv) {
  FLAGS_alsologtostderr = true;
  FLAGS_log_dir = "/tmp";
  google::ParseCommandLineFlags(&argc, &argv, true);
  nglog::InitializeLogging(argv[0]);

  const fs::path certs_dir = ResolvePath(FLAGS_certs_dir);
  CHECK(fs::is_directory(certs_dir))
      << "Not a directory: " << certs_dir.string();

  std::vector<fs::path> files;
  for (const auto &entry : fs::recursive_directory_iterator(certs_dir)) {
    if (entry.is_regular_file()) {
      const std::string entry_path_string = entry.path().string();
      if (entry_path_string.ends_with(".pb.txt") ||
          entry_path_string.ends_with(".pb.txt.gz") ||
          entry_path_string.ends_with(".pb") ||
          entry_path_string.ends_with(".pb.gz")) {
        files.push_back(entry.path());
      }
    }
  }
  std::sort(files.begin(), files.end());

  std::vector<Row> rows;
  for (const fs::path &file : files) {
    const pb::Certificate cert =
        ReadProtoFromFile<pb::Certificate>(file.string());
    if (cert.constrained_tensors_size() == 0) {
      LOG(WARNING) << "Skipping (no constrained_tensors): " << file.string();
      continue;
    }
    const pb::ConstrainedTensor &last =
        cert.constrained_tensors(cert.constrained_tensors_size() - 1);

    Row row;
    row.file = fs::relative(file, certs_dir).string();
    row.problem_name = cert.problem_name();
    row.p = cert.characteristic();
    row.m = cert.extension_degree();
    row.q = IntPow(row.p, row.m);
    row.na = cert.na();
    row.nb = cert.nb();
    row.nc = cert.nc();
    row.lb = last.has_rank_lower_bound()
                 ? std::to_string(last.rank_lower_bound())
                 : "-";
    row.ub = last.has_rank_upper_bound()
                 ? std::to_string(last.rank_upper_bound())
                 : "-";
    rows.push_back(std::move(row));
  }

  std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
    return a.problem_name < b.problem_name;
  });

  const fs::path out_dir = ResolvePath(FLAGS_out_dir);
  CHECK(fs::is_directory(out_dir)) << "Not a directory: " << out_dir.string();
  const fs::path csv_path = out_dir / "cert_summary.csv";
  const fs::path md_path = out_dir / "cert_summary.md";
  WriteCsv(csv_path, rows);
  WriteMarkdown(md_path, rows);

  LOG(INFO) << "Wrote " << rows.size() << " rows to " << csv_path.string()
            << " and " << md_path.string();
  return 0;
}
