// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// dll_react: reactive dynamic lattice liquid simulation on a fully occupied FCC lattice: free
// radical (FRP) or atom transfer radical (ATRP) polymerization of monomers and cross-linkers,
// optionally with solvent and inert precursor chains. See src/dll/README.md.
//
//   dll_react config.ini [key=value ...]
//   dll_react --print-config          (dump the defaults in INI syntax)
//   dll_react --help

#include <exception>
#include <iostream>
#include <ostream>  // std::println overload for std::ostream
#include <print>
#include <string>
#include <vector>

#include "dll/io/config.hpp"
#include "dll/sim/run.hpp"

namespace {

constexpr std::string_view kUsage =
    "usage: dll_react [config.ini] [key=value ...]\n"
    "       dll_react --print-config\n\n"
    "Reactive dynamic lattice liquid simulation (FRP or ATRP) on a fully occupied FCC lattice.\n"
    "Keys given on the command line override the file. Run with --print-config to see every\n"
    "key and its default.\n\n"
    "Please cite the papers listed in README.md under 'Acknowledgements'.\n";

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    dll::Config cfg;
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
        cfg = dll::Config::from_file(a);
      }
    }
    cfg.apply_overrides(overrides);
    (void)dll::run(std::move(cfg), std::cout);
    return 0;
  } catch (const std::exception& e) {
    std::println(std::cerr, "dll_react: {}", e.what());
    return 1;
  }
}
