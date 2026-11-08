// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// bfm_react: reactive Monte Carlo of mono-/divinyl copolymerization on the 3D bond fluctuation
// model (FRP, ATRP). See src/3dbfm/README.md.
//
//   bfm_react config.ini [key=value ...]
//   bfm_react mechanism=frp lattice="32 32 32" n_units=3481 ...
//   bfm_react --print-config          (dump the defaults in INI syntax)
//   bfm_react --help

#include <exception>
#include <iostream>
#include <ostream>  // std::println overload for std::ostream
#include <print>
#include <string>
#include <vector>

#include "bfm/io/config.hpp"
#include "bfm/sim/simulation.hpp"

namespace {

constexpr std::string_view kUsage =
    "usage: bfm_react [config.ini] [key=value ...]\n"
    "       bfm_react --print-config\n\n"
    "Reactive Monte Carlo simulation of radical copolymerization (FRP / ATRP) of\n"
    "mono- and divinyl monomers on the 3D bond fluctuation model. Keys given on the command\n"
    "line override the file. Run with --print-config to see every key and its default.\n\n"
    "Please cite the papers listed in README.md under 'Acknowledgements'.\n";

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    bfm::Config cfg;
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
        cfg = bfm::Config::from_file(a);
      }
    }
    cfg.apply_overrides(overrides);
    cfg.validate();
    auto scheme = bfm::make_scheme(cfg.mechanism, cfg.rates, cfg.competition_rule);
    bfm::Simulation sim(std::move(cfg), std::move(scheme), std::cout);
    sim.run();
    return 0;
  } catch (const std::exception& e) {
    std::println(std::cerr, "bfm_react: {}", e.what());
    return 1;
  }
}
