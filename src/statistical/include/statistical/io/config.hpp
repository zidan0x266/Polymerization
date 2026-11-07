// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration: an INI-style file (`key = value`, optional [sections], `#` comments) plus
// `key=value` command-line overrides. The composition, chemistry and snapshot keys mean the same
// as in the lattice modules, except that the rate constants are rate coefficients here.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "poly/chem/scheme.hpp"

namespace statistical {

struct Config {
  // [system]
  std::size_t n_units{1000000};
  double monovinyl_fraction{0.96};
  double radical_fraction{0.002};

  // [scheme]
  std::string mechanism{"atrp"};
  poly::chem::RateConstants rates{0.0, 0.005, 0.0025, 0.0};
  bool cyclization{true};  ///< false: no reactions between units of the same molecule

  // [run]
  int samples{1};
  int first_sample{1};
  std::uint64_t seed{0};  ///< 0 -> entropy seed, written to run_config.ini
  double stop_conversion{0.96};
  double snapshot_start{0.03};
  double snapshot_step{0.01};
  std::vector<double> extra_snapshots{};
  long log_every{1000000};  ///< reactions between progress reports

  // [output]
  std::filesystem::path output_dir{"."};
  bool write_psf{true};

  /// Load from an INI file; unknown keys are an error so typos cannot silently change runs.
  static Config from_file(const std::filesystem::path& path);
  /// Apply `key=value` overrides (keys as in the file, optionally section-qualified `run.seed`).
  void apply_overrides(const std::vector<std::string>& overrides);
  /// Sanity checks; throws std::invalid_argument with a readable message.
  void validate() const;
  /// The configuration in the same INI syntax (written next to the results).
  [[nodiscard]] std::string to_ini() const;
  /// All conversion levels at which a snapshot is written, ascending.
  [[nodiscard]] std::vector<double> snapshot_levels() const;

  /// Set one key from its textual value; errors are returned, not thrown.
  [[nodiscard]] std::expected<void, std::string> set(std::string_view key, std::string_view value);
};

}  // namespace statistical
