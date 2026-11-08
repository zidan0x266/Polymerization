// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Run configuration: an INI-style file (`key = value`, optional [sections], `#` comments) plus
// `key=value` command-line overrides. The composition, chemistry and snapshot keys mean the
// same as in the 3dbfm module.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "dll/model/initial.hpp"
#include "dll/sim/reaction.hpp"
#include "poly/chem/scheme.hpp"

namespace dll {

struct Config {
  // [system]
  int cells{32};               ///< FCC cubic cells per box edge; 4 cells^3 sites, all occupied
  std::size_t n_units{131072};  ///< reactive units (initiators, monomers, cross-linkers); the rest is solvent
  double monovinyl_fraction{0.96};
  double radical_fraction{0.002};
  std::size_t chain_length{16};
  std::size_t n_chains{0};  ///< inert precursor chains
  InitialConfiguration initial{InitialConfiguration::RandomWalk};
  bool allow_crossing{false};  ///< let chains pass through each other (phantom chains); for comparison only

  // [scheme]
  std::string mechanism{"atrp"};
  poly::chem::RateConstants rates{0.0, 0.005, 0.0025, 0.0};
  poly::chem::CompetitionRule competition_rule{poly::chem::CompetitionRule::Legacy};
  ReactionPartner reaction_partner{ReactionPartner::RandomNeighbor};
  bool cyclization{true};  ///< false: no reactions between units of the same molecule

  // [run]
  int samples{1};
  int first_sample{1};
  std::uint64_t seed{0};  ///< 0 -> entropy seed, written to run_config.ini
  long equilibration_steps{0};
  long steps{1000000000};  ///< maximum number of time steps after the equilibration
  double stop_conversion{0.96};
  double snapshot_start{0.03};
  double snapshot_step{0.01};
  std::vector<double> extra_snapshots{};
  long log_every{1000};

  // [output]
  std::filesystem::path output_dir{"."};
  bool write_pdb{true};
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

}  // namespace dll
