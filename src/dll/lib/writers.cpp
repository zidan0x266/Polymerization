// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/io/writers.hpp"

#include <cctype>
#include <cmath>
#include <format>
#include <fstream>
#include <iterator>
#include <optional>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
std::optional<poly::chem::Species> species_from_code(int code) {
  switch (code) {
    case 1: return poly::chem::Species::Initiator;
    case 2: return poly::chem::Species::Monovinyl;
    case 4: return poly::chem::Species::Divinyl;
    default: return std::nullopt;
  }
}

}  // namespace

void write_dump_frame(std::ostream& stream, const DllSystem& system, long step, bool with_solvent) {
  const std::size_t n = with_solvent ? system.size() : system.solute_beads();
  const double box = kUnit * system.lattice().extent();
  std::string text;
  text.reserve(64 * n + 256);
  auto out = std::back_inserter(text);
  std::format_to(out, "ITEM: TIMESTEP\n{}\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n", step, n);
  for (int axis = 0; axis < 3; ++axis) std::format_to(out, "0 {:.6f}\n", box);
  std::format_to(out, "ITEM: ATOMS id mol type xu yu zu\n");
  for (Bead b = 0; b < n; ++b) {
    const BeadType t = system.type(b);
    const std::size_t mol = t == BeadType::Solvent ? 0 : system.molecule(b) + std::size_t{1};
    const int type = t == BeadType::Chain ? 1 : (t == BeadType::Solvent ? 2 : 3 + system.view(b).legacy_code());
    const Vec3i p = system.unwrapped(b);
    std::format_to(out, "{} {} {} {:.6f} {:.6f} {:.6f}\n", b + 1, mol, type, kUnit * p.x, kUnit * p.y, kUnit * p.z);
  }
  stream.write(text.data(), static_cast<std::streamsize>(text.size()));
}

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

void write_configuration(const std::filesystem::path& path, const DllSystem& system) {
  std::vector<BondPair> reaction_bonds;
  for (Bead b = system.first_unit(); b < system.solute_beads(); ++b) {
    for (const Bead p : system.bonds(b)) {
      if (b < p) reaction_bonds.emplace_back(b, p);
    }
  }
  std::string text;
  auto out = std::back_inserter(text);
  std::format_to(out,
                 "# dll configuration: cells chain_length n_chains n_units n_bonds, then unwrapped lattice points of the\n"
                 "# chain beads, then units (x y z species radicals vinyls; species 1 initiator, 2 monovinyl,\n"
                 "# 4 divinyl), then the bonds formed by reactions (bead numbers); solvent fills the other sites\n"
                 "{} {} {} {} {}\n",
                 system.lattice().cells(), system.chain_length(), system.n_chains(), system.n_units(), reaction_bonds.size());
  for (Bead b = 0; b < system.solute_beads(); ++b) {
    const Vec3i p = system.unwrapped(b);
    if (system.type(b) == BeadType::Chain) {
      std::format_to(out, "{} {} {}\n", p.x, p.y, p.z);
    } else {
      const UnitState s = system.state(b);
      std::format_to(out, "{} {} {} {} {} {}\n", p.x, p.y, p.z, poly::chem::legacy_species_code(s.species), s.radicals, s.vinyls);
    }
  }
  for (const auto& [a, b] : reaction_bonds) std::format_to(out, "{} {}\n", a + 1, b + 1);
  auto file = open_for_write(path);
  file.write(text.data(), static_cast<std::streamsize>(text.size()));
}

DllSystem read_configuration(const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error(std::format("cannot open {}", path.string()));
  std::stringstream in;
  for (std::string line; std::getline(file, line);) {
    if (const auto hash = line.find('#'); hash != std::string::npos) line.erase(hash);
    in << line << '\n';
  }
  const auto fail = [&](std::string_view what) {
    return std::runtime_error(std::format("{}: {}", path.string(), what));
  };
  int cells = 0;
  std::size_t chain_length = 0, n_chains = 0, n_units = 0, n_bonds = 0;
  if (!(in >> cells >> chain_length >> n_chains >> n_units >> n_bonds)) throw fail("bad header");
  if (cells < poly::FccGeometry::kMinCells || cells > poly::FccGeometry::kMaxCells) throw fail("bad number of cells");
  const std::size_t sites = 4 * static_cast<std::size_t>(cells) * static_cast<std::size_t>(cells) * static_cast<std::size_t>(cells);
  if (chain_length * n_chains + n_units > sites || n_bonds > 2 * n_units) throw fail("header does not fit the lattice");
  std::vector<Vec3i> chains(chain_length * n_chains), points(n_units);
  std::vector<UnitState> units(n_units);
  std::vector<BondPair> bonds(n_bonds);
  for (Vec3i& p : chains) {
    if (!(in >> p.x >> p.y >> p.z)) throw fail("truncated chain beads");
  }
  for (std::size_t u = 0; u < n_units; ++u) {
    int species = 0, radicals = 0, vinyls = 0;
    if (!(in >> points[u].x >> points[u].y >> points[u].z >> species >> radicals >> vinyls)) throw fail("truncated units");
    const auto sp = species_from_code(species);
    if (!sp || radicals < 0 || radicals > 2 || vinyls < 0 || vinyls > 2) throw fail(std::format("bad state of unit {}", u + 1));
    units[u] = {*sp, static_cast<std::uint8_t>(radicals), static_cast<std::uint8_t>(vinyls)};
  }
  for (auto& [a, b] : bonds) {
    if (!(in >> a >> b) || a == 0 || b == 0) throw fail("truncated or invalid bonds");
    --a;
    --b;
  }
  try {
    return DllSystem(cells, chain_length, chains, points, units, bonds);
  } catch (const std::invalid_argument& e) {
    throw fail(e.what());
  }
}

}  // namespace dll
