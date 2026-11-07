// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/io/writers.hpp"

#include <cctype>
#include <format>
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

void write_legacy_configuration(const std::filesystem::path& path, const BfmSystem& system) {
  auto out = open_for_write(path);
  const auto& lat = system.lattice();
  std::println(out, "{} {} {} {}", system.size(), lat.lx(), lat.ly(), lat.lz());
  for (Index id = 0; id < system.size(); ++id) {
    const Monomer& m = system.monomer(id);
    const Vec3i a = m.anchor();
    const Vec3i u = system.unwrapped(id);
    std::println(out, "{} {} {} {} {} {}", a.x + 1, a.y + 1, a.z + 1, u.x + 1, u.y + 1, u.z + 1);
    // five partner slots then five bond-vector slots, 1-based, 0 = unused (legacy layout)
    for (int i = 0; i < 5; ++i) std::print(out, "{} ", i < m.nbonds ? m.partners[static_cast<std::size_t>(i)] + 1 : 0);
    for (int i = 0; i < 5; ++i) {
      std::print(out, "{}{}", i < m.nbonds ? m.vectors[static_cast<std::size_t>(i)] + 1 : 0, i < 4 ? ' ' : '\n');
    }
    std::println(out, "{} {} {}", legacy_type_code(m.state, m.nbonds), m.nbonds, legacy_species_code(m.state.species));
  }
}

RestartConfiguration read_legacy_configuration(const std::filesystem::path& path, Vec3i expected_lattice,
                                               std::size_t expected_units) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error(std::format("cannot open {}", path.string()));
  std::size_t n = 0;
  int lx = 0, ly = 0, lz = 0;
  if (!(in >> n >> lx >> ly >> lz)) throw std::runtime_error(std::format("{}: bad header", path.string()));
  if (n != expected_units || Vec3i{lx, ly, lz} != expected_lattice) {
    throw std::runtime_error(std::format("{}: header does not match the configured system size", path.string()));
  }
  RestartConfiguration config{std::vector<Monomer>(n), std::vector<Vec3i>(n)};
  for (std::size_t i = 0; i < n; ++i) {
    Monomer& m = config.units[i];
    int ax, ay, az, ux, uy, uz;
    int partner[5], vec[5];
    int type, nb, species;
    if (!(in >> ax >> ay >> az >> ux >> uy >> uz)) throw std::runtime_error(std::format("{}: truncated at unit {}", path.string(), i));
    for (int& p : partner) in >> p;
    for (int& v : vec) in >> v;
    if (!(in >> type >> nb >> species)) throw std::runtime_error(std::format("{}: truncated at unit {}", path.string(), i));
    if (ax < 1 || ay < 1 || az < 1 || ax > lx || ay > ly || az > lz) {
      throw std::runtime_error(std::format("{}: unit {} lies outside the box", path.string(), i));
    }
    m.set_anchor({ax - 1, ay - 1, az - 1});
    config.unwrapped[i] = {ux - 1, uy - 1, uz - 1};
    if (nb < 0 || nb > kMaxBondsPerUnit) throw std::runtime_error(std::format("{}: bad bond count at unit {}", path.string(), i));
    m.nbonds = static_cast<std::uint8_t>(nb);
    for (int b = 0; b < nb; ++b) {
      if (vec[b] < 1 || vec[b] > kNumBondVectors) throw std::runtime_error(std::format("{}: bad bond index at unit {}", path.string(), i));
      m.partners[static_cast<std::size_t>(b)] = partner[b] - 1;
      m.vectors[static_cast<std::size_t>(b)] = static_cast<BondIndex>(vec[b] - 1);
    }
    // Decode the legacy type code back into (species, radicals, vinyls).
    const Species sp = (species == 1) ? Species::Initiator : (species == 4 ? Species::Divinyl : Species::Monovinyl);
    UnitState s;
    s.species = sp;
    switch (type) {
      case 0: s.radicals = 0; s.vinyls = 1; break;
      case 1: s.radicals = 1; s.vinyls = 0; break;
      case 2: s.radicals = 0; s.vinyls = 0; break;
      case 3: s.radicals = 0; s.vinyls = 2; break;
      case 4: s.radicals = 1; s.vinyls = 1; break;
      case 5: s.radicals = 0; s.vinyls = 1; break;
      case 6: s.radicals = 1; s.vinyls = 0; break;
      case 7: s.radicals = 2; s.vinyls = 0; break;
      case 8: s.radicals = 0; s.vinyls = 0; break;
      case 9: s.radicals = 0; s.vinyls = 0; break;
      default: throw std::runtime_error(std::format("{}: unknown legacy type code {}", path.string(), type));
    }
    m.state = s;
  }
  return config;
}

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

void write_xyz(const std::filesystem::path& path, const BfmSystem& system) {
  auto out = open_for_write(path);
  std::println(out, "{}\nlegacy_type x y z", system.size());
  for (const Monomer& m : system.monomers()) {
    std::println(out, "{:>7} {:>8.3f} {:>8.3f} {:>8.3f}", legacy_type_code(m.state, m.nbonds), static_cast<double>(m.anchor().x + 1),
                 static_cast<double>(m.anchor().y + 1), static_cast<double>(m.anchor().z + 1));
  }
}

}  // namespace bfm
