// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output: PSF snapshots of the network topology in the layout of the lattice modules (legacy
// residue names, each bond listed from both ends), for network analysis. The units have no
// coordinates, so there is no PDB file. The per-sample tables come from src/common
// (poly/chem/tables.hpp).
#pragma once

#include <filesystem>

#include "statistical/model/network.hpp"

namespace statistical {

/// PSF file of all units with all bonds.
void write_psf(const std::filesystem::path& path, const Network& network);

}  // namespace statistical
