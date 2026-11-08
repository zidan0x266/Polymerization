// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/io/writers.hpp"

#include <format>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>

#include "poly/chem/psf.hpp"

namespace bfm {

namespace {

std::ofstream open_for_write(const std::filesystem::path& path, std::ios::openmode mode = std::ios::out) {
  std::ofstream out(path, mode);
  if (!out) throw std::runtime_error(std::format("cannot open {} for writing", path.string()));
  return out;
}

}  // namespace

void write_pdb(const std::filesystem::path& path, const BfmSystem& system) {
  auto out = open_for_write(path);
  const auto& lat = system.lattice();
  std::println(out, "TITLE     DI-VINYL SYSTEM FOR SANOJA");
  std::println(out, "REMARK    SHOW IN SIMULATION BOX");
  std::println(out, "CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P 1", lat.lx() * 1.0, lat.ly() * 1.0, lat.lz() * 1.0, 90.0,
               90.0, 90.0);
  Index id = 0;
  for (const Monomer& m : system.monomers()) {
    const std::string_view res = legacy_residue_name(legacy_type_code(m.state, m.nbonds));
    std::string res_lower(res);
    for (char& c : res_lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    std::println(out, "{:>4}{:>7}  {:>4}{:>4} {:>4}{:>1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}", "ATOM", id + 1, res, res_lower,
                 1, " ", static_cast<double>(m.anchor().x + 1), static_cast<double>(m.anchor().y + 1),
                 static_cast<double>(m.anchor().z + 1), 0.0, 0.0);
    ++id;
  }
  std::println(out, "END");
}

void write_psf(const std::filesystem::path& path, const BfmSystem& system) {
  const auto units = system.monomers();
  poly::chem::write_psf(
      path, units.size(),
      [&](std::size_t i) { return legacy_residue_name(legacy_type_code(units[i].state, units[i].nbonds)); },
      [&](std::size_t i) { return std::span<const Index>(units[i].partners.data(), units[i].nbonds); });
}

}  // namespace bfm
