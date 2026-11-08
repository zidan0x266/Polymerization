// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// rw_pack: packs precursor chains (self-avoiding random walks on an FCC lattice) and single-bead
// cross-linkers into a periodic box and writes a LAMMPS data file for end-linking simulations.
// See src/randomwalk/README.md.
//
//   rw_pack config.ini [key=value ...]
//   rw_pack chain_length=8 n_chains=12416 n_crosslinkers=6208
//   rw_pack --print-config          (dump the defaults in INI syntax)
//   rw_pack --help

#include <exception>
#include <iostream>
#include <ostream>  // std::println overload for std::ostream
#include <print>
#include <string>
#include <vector>

#include "rw/io/config.hpp"
#include "rw/sim/run.hpp"

namespace {

constexpr std::string_view kUsage =
    "usage: rw_pack [config.ini] [key=value ...]\n"
    "       rw_pack --print-config\n\n"
    "Packs precursor chains, grown as self-avoiding random walks on an FCC lattice, and\n"
    "single-bead cross-linkers into a periodic box and writes a LAMMPS data file for\n"
    "end-linking simulations. Keys given on the command line override the file. Run with\n"
    "--print-config to see every key and its default.\n\n"
    "Please cite the papers listed in README.md under 'Acknowledgements'.\n";

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    rw::Config cfg;
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
        cfg = rw::Config::from_file(a);
      }
    }
    cfg.apply_overrides(overrides);
    rw::run(std::move(cfg), std::cout);
    return 0;
  } catch (const std::exception& e) {
    std::println(std::cerr, "rw_pack: {}", e.what());
    return 1;
  }
}
