// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration: an INI-style file (`key = value`, optional [sections], `#` comments) plus
// `key=value` command-line overrides.
#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "dll/model/initial.hpp"

namespace dll {

struct Config {
  // [system]
  int cells{16};  ///< FCC cubic cells per box edge; 4 cells^3 sites, all occupied
  std::size_t chain_length{16};
  std::size_t n_chains{512};  ///< the remaining sites hold solvent
  InitialConfiguration initial{InitialConfiguration::RandomWalk};

  // [run]
  std::uint64_t seed{0};  ///< 0 -> entropy seed, written to run_config.ini
  long steps{10000};
  long log_every{1000};

  // [output]
  std::filesystem::path output_dir{"."};
  long dump_every{0};  ///< write a trajectory frame every this many steps; 0 = none
  bool dump_solvent{false};

  /// Load from an INI file; unknown keys are an error so typos cannot silently change runs.
  static Config from_file(const std::filesystem::path& path);
  /// Apply `key=value` overrides (keys as in the file, optionally section-qualified `run.seed`).
  void apply_overrides(const std::vector<std::string>& overrides);
  /// Sanity checks; throws std::invalid_argument with a readable message.
  void validate() const;
  /// The configuration in the same INI syntax (written next to the results).
  [[nodiscard]] std::string to_ini() const;

  /// Set one key from its textual value; errors are returned, not thrown.
  [[nodiscard]] std::expected<void, std::string> set(std::string_view key, std::string_view value);
};

}  // namespace dll
