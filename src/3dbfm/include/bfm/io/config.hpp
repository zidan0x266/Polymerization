// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration. Replaces the compile-time #defines of the legacy programs with a small
// INI-style file (`key = value`, optional [sections], `#` comments) plus `key=value` command
// line overrides, so one binary serves FRP, ATRP and any lattice size.
#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "bfm/reaction/scheme.hpp"

namespace bfm {

struct Config {
  // [system]
  int lattice_x{100};
  int lattice_y{100};
  int lattice_z{100};
  std::size_t n_units{106250};
  double monovinyl_fraction{0.96};
  double radical_fraction{0.002};

  // [scheme]
  std::string mechanism{"atrp"};
  RateConstants rates{0.0, 0.005, 0.0025, 0.0};
  CompetitionRule competition_rule{CompetitionRule::Legacy};
  /// Legacy behavior: an initiation attempt is made only when a unit happens to sit at the
  /// sampled reaction offset (the initiator branch lived inside the "partner found" block).
  bool initiation_requires_neighbor{true};

  // [run]
  int samples{5};
  int first_sample{1};
  std::uint64_t seed{0};  ///< 0 -> entropy seed, logged for reproducibility
  long mcs_initial_relaxation{1000};
  long mcs_equilibration{10000};
  long max_mcs{990000000};
  double stop_conversion{0.96};
  double snapshot_start{0.03};
  double snapshot_step{0.01};
  std::vector<double> extra_snapshots{};  ///< additional conversion levels (e.g. 0.969)
  long log_every_mcs{1000};
  std::filesystem::path initial_configuration{};  ///< legacy-format restart file; skips placement

  // [output]
  std::filesystem::path output_dir{"."};
  bool write_pdb{true};
  bool write_psf{true};
  bool write_restart{true};

  /// Load from an INI file; unknown keys are an error so typos cannot silently change runs.
  static Config from_file(const std::filesystem::path& path);
  /// Apply `key=value` overrides (keys as in the file, optionally section-qualified `run.seed`).
  void apply_overrides(const std::vector<std::string>& overrides);
  /// Sanity checks; throws std::invalid_argument with a readable message.
  void validate() const;
  /// Human-readable dump in the same INI syntax (written next to the results).
  [[nodiscard]] std::string to_ini() const;

  /// All conversion levels at which a snapshot is written, ascending, deduplicated.
  [[nodiscard]] std::vector<double> snapshot_levels() const;

  /// Set one key from its textual value. Errors are returned, not thrown, so callers can
  /// prefix them with file/line context before reporting.
  [[nodiscard]] std::expected<void, std::string> set(std::string_view key, std::string_view value);
};

}  // namespace bfm
