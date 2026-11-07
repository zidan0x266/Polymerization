// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/io/config.hpp"

#include <format>
#include <sstream>
#include <stdexcept>

#include "poly/ini.hpp"

namespace rw {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_bool;
using poly::ini::parse_number;
using poly::ini::trim;

namespace {

/// Angles in degrees, separated by spaces or commas, e.g. "90, 120".
poly::ini::Parsed<AngleSet> parse_angles(std::string_view key, std::string_view value) {
  std::string text(value);
  for (char& ch : text) {
    if (ch == ',') ch = ' ';
  }
  AngleSet angles;
  std::istringstream in{text};
  std::string token;
  while (in >> token) {
    const auto degrees = parse_number<int>(key, token);
    if (!degrees || !angles.allow(*degrees)) {
      return std::unexpected(std::format("config key '{}': '{}' is not one of 60, 90, 120, 180", key, token));
    }
  }
  if (angles.empty()) return std::unexpected(std::format("config key '{}': no angle given", key));
  return angles;
}

std::string to_string(AngleSet angles) {
  std::string out;
  for (int degrees : AngleSet::kValues) {
    if (angles.allows(degrees)) out += std::format("{}{}", out.empty() ? "" : ", ", degrees);
  }
  return out;
}

}  // namespace

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (key == "chain_length") return assign(chain_length, parse_number<std::size_t>(key, value));
  if (key == "n_chains") return assign(n_chains, parse_number<std::size_t>(key, value));
  if (key == "n_crosslinkers") return assign(n_crosslinkers, parse_number<std::size_t>(key, value));
  if (key == "occupancy") return assign(occupancy, parse_number<double>(key, value));
  if (key == "cells") return assign(cells, parse_number<int>(key, value));
  if (key == "bond_length") return assign(bond_length, parse_number<double>(key, value));
  if (key == "mass") return assign(mass, parse_number<double>(key, value));
  if (key == "bond_angles") return assign(walk.angles, parse_angles(key, value));
  if (key == "site_weighting") {
    const auto w = parse_site_weighting(lower(trim(value)));
    if (!w) return std::unexpected("config key 'site_weighting': expected 'free_neighbors' or 'uniform'");
    walk.weighting = *w;
    return {};
  }
  if (key == "max_attempts") return assign(walk.max_attempts, parse_number<int>(key, value));
  if (key == "seed") return assign(seed, parse_number<std::uint64_t>(key, value));
  if (key == "directory" || key == "output_dir") {
    output_dir = trim(value);
    return {};
  }
  if (key == "data_file") {
    data_file = trim(value);
    return {};
  }
  if (key == "write_xyz") return assign(write_xyz, parse_bool(key, value));
  return std::unexpected(std::format("unknown configuration key '{}'", raw_key));
}

Config Config::from_file(const std::filesystem::path& path) {
  Config cfg;
  poly::ini::read_file(path, [&cfg](std::string_view key, std::string_view value) { return cfg.set(key, value); });
  return cfg;
}

void Config::apply_overrides(const std::vector<std::string>& overrides) {
  poly::ini::apply_overrides(overrides, [this](std::string_view key, std::string_view value) { return set(key, value); });
}

int Config::box_cells() const { return cells > 0 ? cells : cells_for_occupancy(n_chains * chain_length, occupancy); }

void Config::validate() const {
  const auto fail = [](std::string_view msg) { throw std::invalid_argument(std::format("invalid configuration: {}", msg)); };
  if (chain_length < 2) fail("chain_length must be at least 2");
  if (n_chains < 1) fail("n_chains must be at least 1");
  if (cells < 0) fail("cells must be positive, or 0 to derive it from occupancy");
  if (cells == 0 && !(occupancy > 0 && occupancy < 1)) fail("occupancy must be in (0,1)");
  if (!(bond_length > 0)) fail("bond_length must be positive");
  if (!(mass > 0)) fail("mass must be positive");
  if (walk.max_attempts < 1) fail("max_attempts must be at least 1");
  if (data_file.empty()) fail("data_file must not be empty");

  const int n = box_cells();
  if (n < FccLattice::kMinCells || n > FccLattice::kMaxCells) {
    fail(std::format("the box must have {} to {} cells per edge, got {}", FccLattice::kMinCells, FccLattice::kMaxCells, n));
  }
  const double sites = 4.0 * n * n * n;
  if (static_cast<double>(n_chains * chain_length + n_crosslinkers) >= sites) fail("more beads than lattice sites");
}

std::string Config::to_ini() const {
  const auto yesno = [](bool b) { return b ? "true" : "false"; };
  return std::format(
      "[system]\n"
      "chain_length = {}\n"
      "n_chains = {}\n"
      "n_crosslinkers = {}\n"
      "occupancy = {}\n"
      "cells = {}\n"
      "bond_length = {}\n"
      "mass = {}\n\n"
      "[walk]\n"
      "bond_angles = {}\n"
      "site_weighting = {}\n"
      "max_attempts = {}\n\n"
      "[run]\n"
      "seed = {}\n\n"
      "[output]\n"
      "directory = {}\n"
      "data_file = {}\n"
      "write_xyz = {}\n",
      chain_length, n_chains, n_crosslinkers, occupancy, cells, bond_length, mass, to_string(walk.angles),
      rw::to_string(walk.weighting), walk.max_attempts, seed, output_dir.string(), data_file, yesno(write_xyz));
}

}  // namespace rw
