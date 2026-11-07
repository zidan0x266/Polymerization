// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Output of a packing: a LAMMPS data file (atom style molecular, with image flags and atom
// type labels) and an extended XYZ file for visualization.
#pragma once

#include <array>
#include <filesystem>
#include <ostream>
#include <string_view>

#include "rw/model/packing.hpp"

namespace rw {

/// LAMMPS atom type labels, indexed by type - 1.
inline constexpr std::array<std::string_view, kNumBeadTypes> kTypeLabels{"bead", "end", "xlink", "xlink_full"};

/// LAMMPS data file; `title` becomes the first (comment) line. Every type gets mass `mass`.
void write_lammps_data(std::ostream& out, const Packing& packing, double mass, std::string_view title);
void write_lammps_data(const std::filesystem::path& path, const Packing& packing, double mass, std::string_view title);

/// Extended XYZ: atom type as species, wrapped positions, box in the comment line.
void write_xyz(std::ostream& out, const Packing& packing);
void write_xyz(const std::filesystem::path& path, const Packing& packing);

}  // namespace rw
