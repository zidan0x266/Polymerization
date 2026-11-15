// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/io/writers.hpp"

#include <array>
#include <cstddef>
#include <span>

#include "poly/io/pdb.hpp"
#include "poly/io/psf.hpp"

namespace bfm {

void write_pdb(const std::filesystem::path& path, const BfmSystem& system) {
  const auto units = system.monomers();
  const auto& lat = system.lattice();
  poly::io::write_pdb(
      path, "DI-VINYL SYSTEM FOR SANOJA", {lat.lx() * 1.0, lat.ly() * 1.0, lat.lz() * 1.0}, units.size(),
      [&](std::size_t i) { return legacy_residue_name(legacy_type_code(units[i].state, units[i].nbonds)); },
      [&](std::size_t i) {
        const Vec3i a = units[i].anchor();
        return std::array<double, 3>{static_cast<double>(a.x + 1), static_cast<double>(a.y + 1), static_cast<double>(a.z + 1)};
      });
}

void write_psf(const std::filesystem::path& path, const BfmSystem& system) {
  const auto units = system.monomers();
  poly::io::write_psf(
      path, units.size(),
      [&](std::size_t i) { return legacy_residue_name(legacy_type_code(units[i].state, units[i].nbonds)); },
      [&](std::size_t i) { return std::span<const Index>(units[i].partners.data(), units[i].nbonds); });
}

}  // namespace bfm
