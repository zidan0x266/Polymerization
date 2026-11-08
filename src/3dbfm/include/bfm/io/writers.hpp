// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output: PDB/PSF snapshots in the legacy layout (1-based ids and coordinates, legacy residue
// names) for visualization and network analysis. The per-sample tables come from src/common
// (poly/chem/tables.hpp).
#pragma once

#include <filesystem>

#include "bfm/model/system.hpp"
#include "poly/chem/tables.hpp"

namespace bfm {

void write_pdb(const std::filesystem::path& path, const BfmSystem& system);
void write_psf(const std::filesystem::path& path, const BfmSystem& system);

// The per-sample tables are shared with the dll module.
using poly::chem::legacy_conversion_tag;
using poly::chem::SampleTables;

}  // namespace bfm
