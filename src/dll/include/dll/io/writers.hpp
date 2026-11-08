// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output: PDB/PSF snapshots in the layout of the 3dbfm module (legacy residue names, each bond
// listed from both ends), for visualization and network analysis. Lengths are in units of the
// bond length; the solvent is left out. The per-sample tables come from src/common
// (poly/chem/tables.hpp).
#pragma once

#include <filesystem>

#include "dll/model/system.hpp"

namespace dll {

/// PDB file of the chain beads and units at their wrapped positions.
void write_pdb(const std::filesystem::path& path, const DllSystem& system);
/// PSF file of the chain beads and units with all bonds.
void write_psf(const std::filesystem::path& path, const DllSystem& system);

}  // namespace dll
