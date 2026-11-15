// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/io/writers.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <string_view>

#include "poly/io/pdb.hpp"
#include "poly/io/psf.hpp"

namespace dll {

namespace {

const double kUnit = 1.0 / std::sqrt(2.0);  // lattice coordinates to bond lengths

std::string_view residue(const DllSystem& system, Bead b) {
  if (system.type(b) == BeadType::Chain) return "CHN";
  return poly::chem::legacy_residue_name(system.view(b).legacy_code());
}

}  // namespace

void write_pdb(const std::filesystem::path& path, const DllSystem& system) {
  const double box = kUnit * system.lattice().extent();
  poly::io::write_pdb(
      path, "DYNAMIC LATTICE LIQUID SYSTEM", {box, box, box}, system.solute_beads(),
      [&](std::size_t b) { return residue(system, static_cast<Bead>(b)); },
      [&](std::size_t b) {
        const Vec3i p = system.lattice().position(system.site(static_cast<Bead>(b)));
        return std::array<double, 3>{kUnit * p.x, kUnit * p.y, kUnit * p.z};
      });
}

void write_psf(const std::filesystem::path& path, const DllSystem& system) {
  poly::io::write_psf(
      path, system.solute_beads(), [&](std::size_t b) { return residue(system, static_cast<Bead>(b)); },
      [&](std::size_t b) { return system.bonds(static_cast<Bead>(b)); });
}

}  // namespace dll
