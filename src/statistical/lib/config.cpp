// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "statistical/io/config.hpp"

#include <format>
#include <stdexcept>

#include "poly/chem/tables.hpp"
#include "poly/ini.hpp"

namespace statistical {

using poly::ini::assign;
using poly::ini::lower;
using poly::ini::parse_bool;
using poly::ini::parse_double_list;
using poly::ini::parse_number;
using poly::ini::trim;

std::expected<void, std::string> Config::set(std::string_view raw_key, std::string_view value) {
  const std::string key = poly::ini::normalize_key(raw_key);

  if (key == "n_units") return assign(n_units, parse_number<std::size_t>(key, value));
  if (key == "monovinyl_fraction") return assign(monovinyl_fraction, parse_number<double>(key, value));
  if (key == "radical_fraction") return assign(radical_fraction, parse_number<double>(key, value));
  if (key == "mechanism") {
    mechanism = lower(trim(value));
    return {};
  }
  if (key == "k_init") return assign(rates.k_init, parse_number<double>(key, value));
  if (key == "k_p1") return assign(rates.k_p1, parse_number<double>(key, value));
  if (key == "k_px") return assign(rates.k_px, parse_number<double>(key, value));
  if (key == "k_p2") return assign(rates.k_p2, parse_number<double>(key, value));
  if (key == "k_term") return assign(rates.k_term, parse_number<double>(key, value));
  if (key == "cyclization") return assign(cyclization, parse_bool(key, value));
  if (key == "samples") return assign(samples, parse_number<int>(key, value));
  if (key == "first_sample") return assign(first_sample, parse_number<int>(key, value));
  if (key == "seed") return assign(seed, parse_number<std::uint64_t>(key, value));
  if (key == "stop_conversion") return assign(stop_conversion, parse_number<double>(key, value));
  if (key == "snapshot_start") return assign(snapshot_start, parse_number<double>(key, value));
  if (key == "snapshot_step") return assign(snapshot_step, parse_number<double>(key, value));
  if (key == "extra_snapshots") return assign(extra_snapshots, parse_double_list(key, value));
  if (key == "log_every") return assign(log_every, parse_number<long>(key, value));
  if (key == "directory" || key == "output_dir") {
    output_dir = trim(value);
    return {};
  }
  if (key == "write_psf") return assign(write_psf, parse_bool(key, value));
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
  if (n_units < 2 || n_units > 4'000'000'000) fail("n_units must be from 2 to 4e9");
  if (monovinyl_fraction < 0 || monovinyl_fraction > 1) fail("monovinyl_fraction must be in [0,1]");
  if (radical_fraction < 0 || radical_fraction > monovinyl_fraction) fail("radical_fraction must be in [0, monovinyl_fraction]");
  if (static_cast<std::size_t>(static_cast<double>(n_units) * radical_fraction) == 0) fail("radical_fraction yields zero initiators");
  try {
    (void)poly::chem::find_mechanism(mechanism);
  } catch (const std::invalid_argument& e) {
    fail(e.what());
  }
  for (double k : {rates.k_init, rates.k_p1, rates.k_p2, rates.k_term, rates.k_px.value_or(0.0)}) {
    if (!(k >= 0)) fail("rate constants must not be negative");
  }
  if (samples < 1) fail("samples must be >= 1");
  if (first_sample < 0) fail("first_sample must be >= 0");
  if (!(stop_conversion > 0 && stop_conversion <= 1)) fail("stop_conversion must be in (0,1]");
  if (snapshot_step <= 0) fail("snapshot_step must be positive");
  if (snapshot_start < 0 || snapshot_start > 1) fail("snapshot_start must be in [0,1]");
  if (log_every < 1) fail("log_every must be at least 1");
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
      "n_units = {}\n"
      "monovinyl_fraction = {}\n"
      "radical_fraction = {}\n\n"
      "[scheme]\n"
      "mechanism = {}\n"
      "k_init = {}\n"
      "k_p1 = {}\n"
      "{}"
      "k_p2 = {}\n"
      "k_term = {}\n"
      "cyclization = {}\n\n"
      "[run]\n"
      "samples = {}\n"
      "first_sample = {}\n"
      "seed = {}\n"
      "stop_conversion = {}\n"
      "snapshot_start = {}\n"
      "snapshot_step = {}\n"
      "extra_snapshots = {}\n"
      "log_every = {}\n\n"
      "[output]\n"
      "directory = {}\n"
      "write_psf = {}\n",
      n_units, monovinyl_fraction, radical_fraction, mechanism, rates.k_init, rates.k_p1, poly::chem::k_px_ini_line(rates), rates.k_p2,
      rates.k_term, yesno(cyclization), samples, first_sample, seed, stop_conversion, snapshot_start, snapshot_step, extra, log_every,
      output_dir.string(), yesno(write_psf));
}

}  // namespace statistical
