// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration: an INI-style file (`key = value`, optional [sections], `#` comments) plus
// `key=value` command-line overrides. The defaults are the No = 20 system of the fracture paper
// (Macromolecules 2026, 59, 99).
#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "rw/walk/packer.hpp"

namespace rw {

struct Config {
  // [system]
  std::size_t chain_length{18};
  std::size_t n_chains{5680};
  std::size_t n_crosslinkers{2840};
  /// Sizes the box: the smallest number of cells in which the chain beads fill less than this
  /// fraction of the lattice sites (cross-linkers are not counted, as in the original program).
  double occupancy{0.05};
  /// Cells per box edge; 0 derives it from `occupancy`.
  int cells{0};
  double bond_length{1.0};
  double mass{1.0};

  // [walk]
  WalkRules walk{};

  // [run]
  std::uint64_t seed{0};  ///< 0 -> entropy seed, written to run_config.ini

  // [output]
  std::filesystem::path output_dir{"."};
  std::string data_file{"config.lmp"};
  bool write_xyz{true};

  /// Load from an INI file; unknown keys are an error so typos cannot silently change runs.
  static Config from_file(const std::filesystem::path& path);
  /// Apply `key=value` overrides (keys as in the file, optionally section-qualified `run.seed`).
  void apply_overrides(const std::vector<std::string>& overrides);
  /// Sanity checks; throws std::invalid_argument with a readable message.
  void validate() const;
  /// The configuration in the same INI syntax (written next to the results).
  [[nodiscard]] std::string to_ini() const;

  /// Cells per box edge: `cells` if set, else derived from `occupancy`.
  [[nodiscard]] int box_cells() const;

  /// Set one key from its textual value; errors are returned, not thrown.
  [[nodiscard]] std::expected<void, std::string> set(std::string_view key, std::string_view value);
};

}  // namespace rw
