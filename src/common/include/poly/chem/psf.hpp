// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// PSF topology files in the layout of the original 3DBFM code, shared by the modules: one atom
// per unit with its legacy residue name, and every bond listed from both ends.
#pragma once

#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>

namespace poly::chem {

/// Write a PSF file of `n_atoms` atoms: `residue(i)` is the residue name of atom i (zero-based)
/// and `partners(i)` the range of zero-based atoms bonded to it.
template <class Residue, class Partners>
void write_psf(const std::filesystem::path& path, std::size_t n_atoms, Residue&& residue, Partners&& partners) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error(std::format("cannot open {} for writing", path.string()));
  std::size_t total_bonds = 0;
  for (std::size_t i = 0; i < n_atoms; ++i) {
    for ([[maybe_unused]] const auto p : partners(i)) ++total_bonds;
  }
  std::println(out, "{:>8} {:>6}", n_atoms, "!NATOM");
  for (std::size_t i = 0; i < n_atoms; ++i) {
    std::println(out, "{:>8} {:>4} 1   {:>4}  {:>1}    {:>1}      {:>14.6f}{:>14.6f}{:>12}", i + 1, "MAIN", residue(i), "P", "U", 0.01,
                 0.01, 0);
  }
  std::println(out, "");
  std::println(out, "{:>8} {:>13}", total_bonds, "!NBOND: bonds");
  int on_line = 0;
  for (std::size_t i = 0; i < n_atoms; ++i) {
    for (const auto p : partners(i)) {
      std::print(out, "{:>8}{:>8}", i + 1, p + 1);
      if (++on_line == 4) {
        std::println(out, "");
        on_line = 0;
      }
    }
  }
  std::println(out, "END");
}

}  // namespace poly::chem
