// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The generated system: precursor chains followed by single-bead cross-linkers, stored as
// unwrapped lattice points so that bond vectors and periodic image flags are exact.
#pragma once

#include <cstddef>
#include <expected>
#include <numbers>
#include <string>
#include <vector>

#include "rw/lattice/fcc.hpp"

namespace rw {

/// Atom types of the LAMMPS data file. No atom has type 4 at first: the end-linking input
/// assigns it to a cross-linker once all its sites have reacted.
enum class BeadType : int { Interior = 1, End = 2, Crosslinker = 3, ReactedCrosslinker = 4 };
inline constexpr int kNumBeadTypes = 4;

struct Packing {
  int cells{0};
  double bond_length{1.0};
  std::size_t chain_length{0};
  std::size_t n_chains{0};
  std::size_t n_crosslinkers{0};
  /// Unwrapped doubled lattice coordinates: bead v of chain c at c * chain_length + v, then
  /// the cross-linkers.
  std::vector<Vec3i> beads;

  [[nodiscard]] std::size_t size() const noexcept { return beads.size(); }
  [[nodiscard]] std::size_t chain_beads() const noexcept { return n_chains * chain_length; }
  [[nodiscard]] std::size_t n_bonds() const noexcept { return n_chains * (chain_length - 1); }
  [[nodiscard]] bool is_chain_bead(std::size_t i) const noexcept { return i < chain_beads(); }
  /// One-based molecule ID: chain c is c + 1, cross-linker k is n_chains + k + 1.
  [[nodiscard]] std::size_t molecule(std::size_t i) const noexcept {
    return is_chain_bead(i) ? i / chain_length + 1 : i - chain_beads() + n_chains + 1;
  }
  [[nodiscard]] BeadType type(std::size_t i) const noexcept {
    if (!is_chain_bead(i)) return BeadType::Crosslinker;
    const std::size_t v = i % chain_length;
    return (v == 0 || v + 1 == chain_length) ? BeadType::End : BeadType::Interior;
  }
  /// Length of one doubled lattice unit, b/√2.
  [[nodiscard]] double unit() const noexcept { return bond_length / std::numbers::sqrt2; }
  /// Box edge, n a = √2 n b.
  [[nodiscard]] double box_length() const noexcept { return 2 * cells * unit(); }
  /// Period of the doubled coordinates.
  [[nodiscard]] int extent() const noexcept { return 2 * cells; }
};

/// Check the invariants of a packing: one bead per site, every bond a nearest-neighbor step,
/// every bond angle in `angles`, cross-linkers inside the box. Returns the first violation.
[[nodiscard]] std::expected<void, std::string> validate(const Packing& packing, AngleSet angles);

}  // namespace rw
