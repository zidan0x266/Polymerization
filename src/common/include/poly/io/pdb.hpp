// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// PDB coordinate files in the layout of the original 3DBFM code: one atom per unit, named by
// its legacy residue name in upper case, with that name in lower case as the residue name.
#pragma once

#include <array>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <string>
#include <string_view>

#include "poly/io/output_file.hpp"

namespace poly::io {

/// Write a PDB file of `n_atoms` atoms in an orthorhombic box of edges `box`: `residue(i)` is
/// the residue name of atom i (zero-based) and `position(i)` its coordinates as an
/// std::array<double, 3>. `title` fills the TITLE record.
template <class Residue, class Position>
void write_pdb(const std::filesystem::path& path, std::string_view title, const std::array<double, 3>& box, std::size_t n_atoms,
               Residue&& residue, Position&& position) {
  auto out = open_output(path);
  std::println(out, "TITLE     {}", title);
  std::println(out, "REMARK    SHOW IN SIMULATION BOX");
  std::println(out, "CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P 1", box[0], box[1], box[2], 90.0, 90.0, 90.0);
  for (std::size_t i = 0; i < n_atoms; ++i) {
    const std::string_view res = residue(i);
    std::string res_lower(res);
    for (char& c : res_lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    const std::array<double, 3> p = position(i);
    std::println(out, "{:>4}{:>7}  {:>4}{:>4} {:>4}{:>1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}", "ATOM", i + 1, res, res_lower, 1,
                 " ", p[0], p[1], p[2], 0.0, 0.0);
  }
  std::println(out, "END");
}

}  // namespace poly::io
