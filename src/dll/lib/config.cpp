// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/io/config.hpp"

#include <format>
#include <stdexcept>

#include "poly/chem/tables.hpp"
#include "poly/fcc.hpp"
#include "poly/ini.hpp"

namespace dll {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_bool;
using poly::ini::parse_double_list;
using poly::ini::parse_number;
using poly::ini::trim;

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (key == "cells") return assign(cells, parse_number<int>(key, value));
  if (key == "n_units") return assign(n_units, parse_number<std::size_t>(key, value));
  if (key == "monovinyl_fraction") return assign(monovinyl_fraction, parse_number<double>(key, value));
  if (key == "radical_fraction") return assign(radical_fraction, parse_number<double>(key, value));
  if (key == "chain_length") return assign(chain_length, parse_number<std::size_t>(key, value));
  if (key == "n_chains") return assign(n_chains, parse_number<std::size_t>(key, value));
  if (key == "initial") {
    const auto c = parse_initial_configuration(lower(trim(value)));
    if (!c) return std::unexpected("config key 'initial': expected 'random_walk' or 'aligned'");
    initial = *c;
    return {};
  }
  if (key == "allow_crossing") return assign(allow_crossing, parse_bool(key, value));
  if (key == "mechanism") {
    mechanism = lower(trim(value));
    return {};
  }
  if (key == "k_init") return assign(rates.k_init, parse_number<double>(key, value));
  if (key == "k_p1") return assign(rates.k_p1, parse_number<double>(key, value));
  if (key == "k_p2") return assign(rates.k_p2, parse_number<double>(key, value));
  if (key == "k_term") return assign(rates.k_term, parse_number<double>(key, value));
  if (key == "competition_rule") {
    const auto r = poly::chem::parse_competition_rule(trim(value));
    if (!r) return std::unexpected("config key 'competition_rule': expected 'legacy' or 'partition'");
    competition_rule = *r;
    return {};
  }
  if (key == "reaction_partner") {
    const auto p = parse_reaction_partner(lower(trim(value)));
    if (!p) return std::unexpected("config key 'reaction_partner': expected 'random_neighbor' or 'all_neighbors'");
    reaction_partner = *p;
    return {};
  }
  if (key == "samples") return assign(samples, parse_number<int>(key, value));
  if (key == "first_sample") return assign(first_sample, parse_number<int>(key, value));
  if (key == "seed") return assign(seed, parse_number<std::uint64_t>(key, value));
  if (key == "equilibration_steps") return assign(equilibration_steps, parse_number<long>(key, value));
  if (key == "steps") return assign(steps, parse_number<long>(key, value));
  if (key == "stop_conversion") return assign(stop_conversion, parse_number<double>(key, value));
  if (key == "snapshot_start") return assign(snapshot_start, parse_number<double>(key, value));
  if (key == "snapshot_step") return assign(snapshot_step, parse_number<double>(key, value));
  if (key == "extra_snapshots") return assign(extra_snapshots, parse_double_list(key, value));
  if (key == "log_every") return assign(log_every, parse_number<long>(key, value));
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
  if (key == "dump_every") return assign(dump_every, parse_number<long>(key, value));
  if (key == "dump_solvent") return assign(dump_solvent, parse_bool(key, value));
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
  if (cells < poly::FccGeometry::kMinCells || cells > poly::FccGeometry::kMaxCells) {
    fail(std::format("cells must be from {} to {}", poly::FccGeometry::kMinCells, poly::FccGeometry::kMaxCells));
  }
  if (n_chains > 0 && chain_length < 2) fail("chain_length must be at least 2");
  const std::size_t sites = 4 * static_cast<std::size_t>(cells) * static_cast<std::size_t>(cells) * static_cast<std::size_t>(cells);
  if (initial_configuration.empty() && n_units + n_chains * chain_length > sites) {
    fail(std::format("n_units and the chain beads exceed the {} sites of {} cells", sites, cells));
  }
  if (monovinyl_fraction < 0 || monovinyl_fraction > 1) fail("monovinyl_fraction must be in [0,1]");
  if (radical_fraction < 0 || radical_fraction > monovinyl_fraction) fail("radical_fraction must be in [0, monovinyl_fraction]");
  if (initial_configuration.empty() && n_units > 0 && static_cast<std::size_t>(static_cast<double>(n_units) * radical_fraction) == 0) {
    fail("radical_fraction yields zero initiators");
  }
  for (double k : {rates.k_init, rates.k_p1, rates.k_p2, rates.k_term}) {
    if (k < 0 || k > 1) fail("rate constants are probabilities in [0,1]");
  }
  if (samples < 1) fail("samples must be >= 1");
  if (first_sample < 0) fail("first_sample must be >= 0");
  if (equilibration_steps < 0 || steps < 0) fail("step counts must not be negative");
  if (!(stop_conversion > 0 && stop_conversion <= 1)) fail("stop_conversion must be in (0,1]");
  if (snapshot_step <= 0) fail("snapshot_step must be positive");
  if (snapshot_start < 0 || snapshot_start > 1) fail("snapshot_start must be in [0,1]");
  if (log_every < 1) fail("log_every must be at least 1");
  if (dump_every < 0) fail("dump_every must not be negative");
}

std::vector<double> Config::snapshot_levels() const {
  return poly::chem::snapshot_levels(snapshot_start, snapshot_step, extra_snapshots, stop_conversion);
}

std::string Config::to_ini() const {
  std::string extra;
  for (std::size_t i = 0; i < extra_snapshots.size(); ++i) extra += std::format("{}{}", i ? ", " : "", extra_snapshots[i]);
  const auto yesno = [](bool b) { return b ? "true" : "false"; };
  return std::format(
      "[system]\n"
      "cells = {}\n"
      "n_units = {}\n"
      "monovinyl_fraction = {}\n"
      "radical_fraction = {}\n"
      "chain_length = {}\n"
      "n_chains = {}\n"
      "initial = {}\n"
      "allow_crossing = {}\n\n"
      "[scheme]\n"
      "mechanism = {}\n"
      "k_init = {}\n"
      "k_p1 = {}\n"
      "k_p2 = {}\n"
      "k_term = {}\n"
      "competition_rule = {}\n"
      "reaction_partner = {}\n\n"
      "[run]\n"
      "samples = {}\n"
      "first_sample = {}\n"
      "seed = {}\n"
      "equilibration_steps = {}\n"
      "steps = {}\n"
      "stop_conversion = {}\n"
      "snapshot_start = {}\n"
      "snapshot_step = {}\n"
      "extra_snapshots = {}\n"
      "log_every = {}\n"
      "initial_configuration = {}\n\n"
      "[output]\n"
      "directory = {}\n"
      "write_pdb = {}\n"
      "write_psf = {}\n"
      "write_restart = {}\n"
      "dump_every = {}\n"
      "dump_solvent = {}\n",
      cells, n_units, monovinyl_fraction, radical_fraction, chain_length, n_chains, to_string(initial), yesno(allow_crossing), mechanism, rates.k_init,
      rates.k_p1, rates.k_p2, rates.k_term, poly::chem::to_string(competition_rule), to_string(reaction_partner), samples,
      first_sample, seed, equilibration_steps, steps, stop_conversion, snapshot_start, snapshot_step, extra, log_every,
      initial_configuration.string(), output_dir.string(), yesno(write_pdb), yesno(write_psf), yesno(write_restart), dump_every,
      yesno(dump_solvent));
}

}  // namespace dll
