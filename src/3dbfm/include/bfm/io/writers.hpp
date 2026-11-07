// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output: PDB/PSF snapshots in the legacy layout (1-based ids and coordinates, legacy residue
// names) and the legacy `ini.dat` restart format for interoperability with existing analysis
// scripts. The per-sample tables come from src/common (poly/chem/tables.hpp).
#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <vector>

#include "bfm/model/system.hpp"
#include "poly/chem/tables.hpp"

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

// The per-sample tables are shared with the dll module.
using poly::chem::legacy_conversion_tag;
using poly::chem::SampleTables;

}  // namespace bfm
