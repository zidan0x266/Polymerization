// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "poly/chem/tables.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>

namespace poly::chem {

std::string legacy_conversion_tag(double conversion) { return std::format("{:f}", conversion); }

SampleTables::SampleTables(const std::filesystem::path& dir, int sample, std::span<const std::string> category_labels,
                           std::string_view time_unit) {
  const auto open = [&](std::string_view stem, std::string_view header) {
    const auto path = dir / std::format("{}_{}.dat", stem, sample);
    const bool fresh = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
    std::ofstream f(path, std::ios::app);
    if (!f) throw std::runtime_error(std::format("cannot open {} for writing", path.string()));
    if (fresh) std::println(f, "{}", header);
    return f;
  };
  structural_ = open("structural", std::format("# {} conversion s0_monomer s1_radical s2_dead s3_crosslinker s4 s5 s6 s7 s8 s9_initiator",
                                                time_unit));
  propagation_ = open("propagation", "# conversion monomer_monomer polymer_monomer polymer_polymer intramolecular");
  std::string header = std::format("# {} conversion", time_unit);
  for (const std::string& label : category_labels) header += std::format(" {}", label);
  reactype_ = open("reactype", header);
}

void SampleTables::append(long time, double conversion, const StructuralCounts& structure, const ConnectivityTally& tally,
                          std::span<const std::uint64_t> category_counts) {
  std::print(structural_, "{} {:f}", time, conversion);
  for (auto c : structure) std::print(structural_, " {}", c);
  std::println(structural_, "");
  std::println(propagation_, "{:f} {} {} {} {}", conversion, tally.monomer_monomer, tally.polymer_monomer, tally.polymer_polymer,
               tally.intramolecular);
  std::print(reactype_, "{} {:f}", time, conversion);
  for (auto c : category_counts) std::print(reactype_, " {}", c);
  std::println(reactype_, "");
  for (std::ofstream* f : {&structural_, &propagation_, &reactype_}) f->flush();
}

std::vector<double> snapshot_levels(double start, double step, std::span<const double> extra, double stop) {
  std::vector<double> levels;
  // integer counter instead of repeated addition: no floating-point drift (the legacy
  // 0.91 + 5*0.01 > 0.96 accident decided whether a 0.969 snapshot was written)
  for (long i = 0;; ++i) {
    const double level = start + static_cast<double>(i) * step;
    if (level > stop + 1e-9) break;
    levels.push_back(level);
  }
  for (double x : extra) {
    if (x > 0 && x <= stop + 1e-9) levels.push_back(x);
  }
  levels.push_back(stop);
  std::ranges::sort(levels);
  const auto dup = std::ranges::unique(levels, [](double a, double b) { return std::fabs(a - b) < 1e-9; });
  levels.erase(dup.begin(), dup.end());
  return levels;
}

}  // namespace poly::chem
