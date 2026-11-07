// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// stat_react: coordinate-free kinetic Monte Carlo of free radical (FRP) and atom transfer radical
// (ATRP) polymerization of monomers and cross-linkers. See src/statistical/README.md.
//
//   stat_react config.ini [key=value ...]
//   stat_react --print-config          (dump the defaults in INI syntax)
//   stat_react --help

#include <exception>
#include <iostream>
#include <ostream>  // std::println overload for std::ostream
#include <print>
#include <string>
#include <vector>

#include "statistical/io/config.hpp"
#include "statistical/sim/run.hpp"

namespace {

constexpr std::string_view kUsage =
    "usage: stat_react [config.ini] [key=value ...]\n"
    "       stat_react --print-config\n\n"
    "Coordinate-free kinetic Monte Carlo of radical polymerization (FRP or ATRP). Keys given on the\n"
    "command line override the file. Run with --print-config to see every key and its default.\n";

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    statistical::Config cfg;
    std::vector<std::string> overrides;
    for (const std::string& a : args) {
      if (a == "-h" || a == "--help") {
        std::print("{}", kUsage);
        return 0;
      }
      if (a == "--print-config") {
        std::print("{}", cfg.to_ini());
        return 0;
      }
      if (a.contains('=')) {
        overrides.push_back(a);
      } else {
        cfg = statistical::Config::from_file(a);
      }
    }
    cfg.apply_overrides(overrides);
    (void)statistical::run(std::move(cfg), std::cout);
    return 0;
  } catch (const std::exception& e) {
    std::println(std::cerr, "stat_react: {}", e.what());
    return 1;
  }
}
