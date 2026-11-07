// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output: per-sample tables (same columns as the legacy structural_/propagation_/reactype_
// files, now with a '#' header line), PDB/PSF snapshots in the legacy layout (1-based ids and
// coordinates, legacy residue names) and the legacy `ini.dat` restart format for
// interoperability with existing analysis scripts.
#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <vector>

#include "bfm/model/system.hpp"

namespace bfm {

/// Units read from a restart file, with their unwrapped anchors (for BfmSystem::load).
struct RestartConfiguration {
  std::vector<Monomer> units;
  std::vector<Vec3i> unwrapped;
};

/// Legacy restart format ("ini.dat"): header `ntot lx ly lz`, then three lines per unit.
void write_legacy_configuration(const std::filesystem::path& path, const BfmSystem& system);
[[nodiscard]] RestartConfiguration read_legacy_configuration(const std::filesystem::path& path, Vec3i expected_lattice,
                                                             std::size_t expected_units);

void write_pdb(const std::filesystem::path& path, const BfmSystem& system);
void write_psf(const std::filesystem::path& path, const BfmSystem& system);
void write_xyz(const std::filesystem::path& path, const BfmSystem& system);

/// Append-mode tables for one sample. Files are created with a header on first use.
class SampleTables {
 public:
  SampleTables(const std::filesystem::path& dir, int sample, std::span<const std::string> category_labels);

  void append(long mcs, double conversion, const BfmSystem& system, std::span<const std::uint64_t> category_counts);

 private:
  std::ofstream structural_;
  std::ofstream propagation_;
  std::ofstream reactype_;
};

/// Format a conversion level the way the legacy file names did ("%lf" -> 6 decimals).
[[nodiscard]] std::string legacy_conversion_tag(double conversion);

}  // namespace bfm
