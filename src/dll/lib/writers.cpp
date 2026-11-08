// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/io/writers.hpp"

#include <cctype>
#include <cmath>
#include <format>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>

#include "poly/chem/psf.hpp"

namespace dll {

namespace {

const double kUnit = 1.0 / std::sqrt(2.0);  // lattice coordinates to bond lengths

std::ofstream open_for_write(const std::filesystem::path& path) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error(std::format("cannot open {} for writing", path.string()));
  return out;
}

std::string_view residue(const DllSystem& system, Bead b) {
  if (system.type(b) == BeadType::Chain) return "CHN";
  return poly::chem::legacy_residue_name(system.view(b).legacy_code());
}

/// Species codes of the configuration file: the legacy `monreactive` column.

}  // namespace

void write_pdb(const std::filesystem::path& path, const DllSystem& system) {
  auto out = open_for_write(path);
  const double box = kUnit * system.lattice().extent();
  std::println(out, "TITLE     DYNAMIC LATTICE LIQUID SYSTEM");
  std::println(out, "REMARK    SHOW IN SIMULATION BOX");
  std::println(out, "CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P 1", box, box, box, 90.0, 90.0, 90.0);
  for (Bead b = 0; b < system.solute_beads(); ++b) {
    const std::string_view res = residue(system, b);
    std::string res_lower(res);
    for (char& c : res_lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    const Vec3i p = system.lattice().position(system.site(b));
    std::println(out, "{:>4}{:>7}  {:>4}{:>4} {:>4}{:>1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}", "ATOM", b + 1, res, res_lower, 1,
                 " ", kUnit * p.x, kUnit * p.y, kUnit * p.z, 0.0, 0.0);
  }
  std::println(out, "END");
}

void write_psf(const std::filesystem::path& path, const DllSystem& system) {
  poly::chem::write_psf(
      path, system.solute_beads(), [&](std::size_t b) { return residue(system, static_cast<Bead>(b)); },
      [&](std::size_t b) { return system.bonds(static_cast<Bead>(b)); });
}

}  // namespace dll
