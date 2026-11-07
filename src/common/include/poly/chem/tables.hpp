// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Per-sample tables written at every snapshot, the same in every lattice module (the legacy
// structural_/propagation_/reactype_ files, now with a '#' header line), and the conversion
// levels at which snapshots are taken.
#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "poly/chem/progress.hpp"
#include "poly/chem/unit.hpp"

namespace poly::chem {

/// Number of units carrying each legacy structural-unit code 0..9 (the legacy s_0..s_9).
using StructuralCounts = std::array<std::int64_t, kNumLegacyTypeCodes>;

/// Append-mode tables for one sample. Files are created with a header on first use; the first
/// column is the time in `time_unit` (MCS in 3dbfm, time steps in dll, a continuous time in the
/// statistical model).
class SampleTables {
 public:
  SampleTables(const std::filesystem::path& dir, int sample, std::span<const std::string> category_labels,
               std::string_view time_unit = "mcs");

  void append(long time, double conversion, const StructuralCounts& structure, const ConnectivityTally& tally,
              std::span<const std::uint64_t> category_counts);
  void append(double time, double conversion, const StructuralCounts& structure, const ConnectivityTally& tally,
              std::span<const std::uint64_t> category_counts);

 private:
  void append_rows(std::string_view time, double conversion, const StructuralCounts& structure, const ConnectivityTally& tally,
                   std::span<const std::uint64_t> category_counts);

  std::ofstream structural_;
  std::ofstream propagation_;
  std::ofstream reactype_;
};

/// Format a conversion level the way the legacy file names did ("%lf" -> 6 decimals).
[[nodiscard]] std::string legacy_conversion_tag(double conversion);

/// All conversion levels at which a snapshot is written: start, start + step, ... up to `stop`,
/// plus the `extra` levels and `stop` itself, ascending and deduplicated.
[[nodiscard]] std::vector<double> snapshot_levels(double start, double step, std::span<const double> extra, double stop);

}  // namespace poly::chem
