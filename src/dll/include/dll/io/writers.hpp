// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output. Trajectories in the LAMMPS dump format (readable by OVITO, VMD and most analysis
// tools), PDB/PSF snapshots in the layout of the 3dbfm module (legacy residue names, each bond
// listed from both ends), and a configuration file from which a run can be continued. Lengths
// are in units of the bond length; the solvent is left out of the snapshots.
#pragma once

#include <filesystem>
#include <ostream>

#include "dll/model/system.hpp"

namespace dll {

/// Append one frame with the chain beads and units, and the solvent if `with_solvent`. Columns:
/// id, mol (a bead of the molecule, 0 for solvent), type (1 chain, 2 solvent, 3 + legacy code
/// for units) and the unwrapped position.
void write_dump_frame(std::ostream& out, const DllSystem& system, long step, bool with_solvent);

/// PDB file of the chain beads and units at their wrapped positions.
void write_pdb(const std::filesystem::path& path, const DllSystem& system);
/// PSF file of the chain beads and units with all bonds.
void write_psf(const std::filesystem::path& path, const DllSystem& system);

/// Configuration file: chains, units with their states, and the bonds formed by reactions.
void write_configuration(const std::filesystem::path& path, const DllSystem& system);
[[nodiscard]] DllSystem read_configuration(const std::filesystem::path& path);

}  // namespace dll
