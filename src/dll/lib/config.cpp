// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/io/config.hpp"

#include <format>
#include <stdexcept>

#include "poly/fcc.hpp"
#include "poly/ini.hpp"

namespace dll {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_bool;
using poly::ini::parse_number;
using poly::ini::trim;

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (key == "cells") return assign(cells, parse_number<int>(key, value));
  if (key == "chain_length") return assign(chain_length, parse_number<std::size_t>(key, value));
  if (key == "n_chains") return assign(n_chains, parse_number<std::size_t>(key, value));
  if (key == "initial") {
    const auto c = parse_initial_configuration(lower(trim(value)));
    if (!c) return std::unexpected("config key 'initial': expected 'random_walk' or 'aligned'");
    initial = *c;
    return {};
  }
  if (key == "seed") return assign(seed, parse_number<std::uint64_t>(key, value));
  if (key == "steps") return assign(steps, parse_number<long>(key, value));
  if (key == "log_every") return assign(log_every, parse_number<long>(key, value));
  if (key == "directory" || key == "output_dir") {
    output_dir = trim(value);
    return {};
  }
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
  const double sites = 4.0 * cells * cells * cells;
  if (static_cast<double>(n_chains * chain_length) > sites) fail("more chain beads than lattice sites");
  if (steps < 0) fail("steps must not be negative");
  if (log_every < 1) fail("log_every must be at least 1");
  if (dump_every < 0) fail("dump_every must not be negative");
}

std::string Config::to_ini() const {
  const auto yesno = [](bool b) { return b ? "true" : "false"; };
  return std::format(
      "[system]\n"
      "cells = {}\n"
      "chain_length = {}\n"
      "n_chains = {}\n"
      "initial = {}\n\n"
      "[run]\n"
      "seed = {}\n"
      "steps = {}\n"
      "log_every = {}\n\n"
      "[output]\n"
      "directory = {}\n"
      "dump_every = {}\n"
      "dump_solvent = {}\n",
      cells, chain_length, n_chains, to_string(initial), seed, steps, log_every, output_dir.string(), dump_every,
      yesno(dump_solvent));
}

}  // namespace dll
