// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/io/config.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <sstream>
#include <stdexcept>

#include "poly/ini.hpp"

namespace bfm {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_bool;
using poly::ini::parse_double_list;
using poly::ini::parse_number;
using poly::ini::trim;

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (key == "lattice") {
    std::istringstream in{trim(value)};
    if (!(in >> lattice_x >> lattice_y >> lattice_z)) return std::unexpected("config key 'lattice': expected three integers");
    return {};
  }
  if (key == "lattice_x") return assign(lattice_x, parse_number<int>(key, value));
  if (key == "lattice_y") return assign(lattice_y, parse_number<int>(key, value));
  if (key == "lattice_z") return assign(lattice_z, parse_number<int>(key, value));
  if (key == "n_units") return assign(n_units, parse_number<std::size_t>(key, value));
  if (key == "monovinyl_fraction") return assign(monovinyl_fraction, parse_number<double>(key, value));
  if (key == "radical_fraction") return assign(radical_fraction, parse_number<double>(key, value));
  if (key == "mechanism" || key == "type" || key == "scheme") {
    mechanism = lower(trim(value));
    return {};
  }
  if (key == "k_init") return assign(rates.k_init, parse_number<double>(key, value));
  if (key == "k_p1") return assign(rates.k_p1, parse_number<double>(key, value));
  if (key == "k_p2") return assign(rates.k_p2, parse_number<double>(key, value));
  if (key == "k_term") return assign(rates.k_term, parse_number<double>(key, value));
  if (key == "competition_rule") {
    const auto r = parse_competition_rule(trim(value));
    if (!r) return std::unexpected("config key 'competition_rule': expected 'legacy' or 'partition'");
    competition_rule = *r;
    return {};
  }
  if (key == "initiation_requires_neighbor" || key == "initiation_requires_neighbour") {  // British spelling accepted
    return assign(initiation_requires_neighbor, parse_bool(key, value));
  }
  if (key == "samples") return assign(samples, parse_number<int>(key, value));
  if (key == "first_sample") return assign(first_sample, parse_number<int>(key, value));
  if (key == "seed") return assign(seed, parse_number<std::uint64_t>(key, value));
  if (key == "mcs_initial_relaxation") return assign(mcs_initial_relaxation, parse_number<long>(key, value));
  if (key == "mcs_equilibration") return assign(mcs_equilibration, parse_number<long>(key, value));
  if (key == "max_mcs") return assign(max_mcs, parse_number<long>(key, value));
  if (key == "stop_conversion") return assign(stop_conversion, parse_number<double>(key, value));
  if (key == "snapshot_start") return assign(snapshot_start, parse_number<double>(key, value));
  if (key == "snapshot_step") return assign(snapshot_step, parse_number<double>(key, value));
  if (key == "extra_snapshots") return assign(extra_snapshots, parse_double_list(key, value));
  if (key == "log_every_mcs") return assign(log_every_mcs, parse_number<long>(key, value));
  if (key == "initial_configuration") {
    initial_configuration = trim(value);
    return {};
  }
  if (key == "directory" || key == "output_dir") {
    output_dir = trim(value);
    return {};
  }
  if (key == "write_pdb") return assign(write_pdb, parse_bool(key, value));
  if (key == "write_psf") return assign(write_psf, parse_bool(key, value));
  if (key == "write_restart") return assign(write_restart, parse_bool(key, value));
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

void Config::validate() const {
  const auto fail = [](std::string_view msg) { throw std::invalid_argument(std::format("invalid configuration: {}", msg)); };
  for (int L : {lattice_x, lattice_y, lattice_z}) {
    if (L < 4 || L % 2 != 0) fail("lattice dimensions must be even and >= 4");
  }
  const std::size_t capacity = static_cast<std::size_t>(lattice_x / 2) * static_cast<std::size_t>(lattice_y / 2) *
                               static_cast<std::size_t>(lattice_z / 2);
  if (n_units == 0) fail("n_units must be positive");
  if (n_units > capacity) fail("n_units exceeds the L^3/8 cubes that fit on the even sub-lattice");
  if (monovinyl_fraction < 0 || monovinyl_fraction > 1) fail("monovinyl_fraction must be in [0,1]");
  if (radical_fraction < 0 || radical_fraction > monovinyl_fraction) fail("radical_fraction must be in [0, monovinyl_fraction]");
  if (static_cast<std::size_t>(static_cast<double>(n_units) * radical_fraction) == 0) fail("radical_fraction yields zero initiators");
  if (samples < 1) fail("samples must be >= 1");
  if (first_sample < 0) fail("first_sample must be >= 0");
  if (mcs_initial_relaxation < 0 || mcs_equilibration < 0 || max_mcs < 1) fail("MCS counts must be non-negative (max_mcs >= 1)");
  if (!(stop_conversion > 0 && stop_conversion <= 1)) fail("stop_conversion must be in (0,1]");
  if (snapshot_step <= 0) fail("snapshot_step must be positive");
  if (snapshot_start < 0 || snapshot_start > 1) fail("snapshot_start must be in [0,1]");
  if (log_every_mcs < 1) fail("log_every_mcs must be >= 1");
  for (double k : {rates.k_init, rates.k_p1, rates.k_p2, rates.k_term}) {
    if (k < 0 || k > 1) fail("rate constants are per-attempt probabilities in [0,1]");
  }
}

std::vector<double> Config::snapshot_levels() const {
  std::vector<double> levels;
  // integer counter instead of repeated addition: no floating-point drift (the legacy
  // 0.91 + 5*0.01 > 0.96 accident decided whether a 0.969 snapshot was written)
  for (long i = 0;; ++i) {
    const double level = snapshot_start + static_cast<double>(i) * snapshot_step;
    if (level > stop_conversion + 1e-9) break;
    levels.push_back(level);
  }
  for (double x : extra_snapshots) {
    if (x > 0 && x <= stop_conversion + 1e-9) levels.push_back(x);
  }
  levels.push_back(stop_conversion);
  std::ranges::sort(levels);
  const auto dup = std::ranges::unique(levels, [](double a, double b) { return std::fabs(a - b) < 1e-9; });
  levels.erase(dup.begin(), dup.end());
  return levels;
}

std::string Config::to_ini() const {
  std::string extra;
  for (std::size_t i = 0; i < extra_snapshots.size(); ++i) extra += std::format("{}{}", i ? ", " : "", extra_snapshots[i]);
  const auto yesno = [](bool b) { return b ? "true" : "false"; };
  return std::format(
      "[system]\n"
      "lattice = {} {} {}\n"
      "n_units = {}\n"
      "monovinyl_fraction = {}\n"
      "radical_fraction = {}\n\n"
      "[scheme]\n"
      "mechanism = {}\n"
      "k_init = {}\n"
      "k_p1 = {}\n"
      "k_p2 = {}\n"
      "k_term = {}\n"
      "competition_rule = {}\n"
      "initiation_requires_neighbor = {}\n\n"
      "[run]\n"
      "samples = {}\n"
      "first_sample = {}\n"
      "seed = {}\n"
      "mcs_initial_relaxation = {}\n"
      "mcs_equilibration = {}\n"
      "max_mcs = {}\n"
      "stop_conversion = {}\n"
      "snapshot_start = {}\n"
      "snapshot_step = {}\n"
      "extra_snapshots = {}\n"
      "log_every_mcs = {}\n"
      "initial_configuration = {}\n\n"
      "[output]\n"
      "directory = {}\n"
      "write_pdb = {}\n"
      "write_psf = {}\n"
      "write_restart = {}\n",
      lattice_x, lattice_y, lattice_z, n_units, monovinyl_fraction, radical_fraction, mechanism, rates.k_init, rates.k_p1,
      rates.k_p2, rates.k_term, to_string(competition_rule), yesno(initiation_requires_neighbor), samples, first_sample, seed,
      mcs_initial_relaxation, mcs_equilibration, max_mcs, stop_conversion, snapshot_start, snapshot_step, extra, log_every_mcs,
      initial_configuration.string(), output_dir.string(), yesno(write_pdb), yesno(write_psf), yesno(write_restart));
}

}  // namespace bfm
