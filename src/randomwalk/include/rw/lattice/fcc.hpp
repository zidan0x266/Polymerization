// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Periodic face-centered cubic (FCC) lattice with site occupancy, on the shared geometry of
// src/common (doubled integer coordinates, twelve nearest neighbors), and the bond-angle rule
// of the random walks.
#pragma once

#include <array>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include "poly/fcc.hpp"

namespace rw {

using poly::Site;
using poly::Vec3i;

/// Set of bond angles a chain may form. Two consecutive nearest-neighbor steps p and q meet at
/// an angle of 180, 120, 90, 60 or 0 degrees for p.q = 2, 1, 0, -1, -2; the 0-degree reversal
/// would put two beads on one site and is never allowed.
class AngleSet {
 public:
  static constexpr std::array<int, 4> kValues{60, 90, 120, 180};

  constexpr AngleSet() = default;
  constexpr AngleSet(std::initializer_list<int> degrees) {
    for (int d : degrees) allow(d);
  }

  /// Add an angle in degrees; returns false if it is not one of kValues.
  constexpr bool allow(int degrees) noexcept {
    for (int k = 0; k < static_cast<int>(kValues.size()); ++k) {
      if (kValues[k] == degrees) {
        mask_ |= static_cast<std::uint8_t>(1u << (k + 1));
        return true;
      }
    }
    return false;
  }
  [[nodiscard]] constexpr bool allows(int degrees) const noexcept {
    for (int k = 0; k < static_cast<int>(kValues.size()); ++k) {
      if (kValues[k] == degrees) return (mask_ >> (k + 1)) & 1u;
    }
    return false;
  }
  /// Whether step `next` may follow step `previous`.
  [[nodiscard]] constexpr bool allows_steps(Vec3i previous, Vec3i next) const noexcept {
    return (mask_ >> (previous.dot(next) + 2)) & 1u;
  }
  [[nodiscard]] constexpr bool empty() const noexcept { return mask_ == 0; }
  constexpr bool operator==(const AngleSet&) const noexcept = default;

 private:
  std::uint8_t mask_{0};  ///< bit p.q + 2
};

/// FCC lattice whose sites are either empty or occupied; each site also counts its occupied
/// nearest neighbors.
class FccLattice : public poly::FccGeometry {
 public:
  /// A lattice of `cells` cubic cells per edge (kMinCells to kMaxCells), all sites empty.
  explicit FccLattice(int cells) : FccGeometry(cells), state_(num_sites(), 0) {}

  [[nodiscard]] bool occupied(Site s) const noexcept { return (state_[s] & kOccupiedBit) != 0; }
  /// Number of occupied nearest neighbors of `s` (0..12).
  [[nodiscard]] int occupied_neighbors(Site s) const noexcept { return state_[s] & kCountMask; }
  [[nodiscard]] std::size_t occupied_count() const noexcept { return occupied_count_; }

  /// Mark `s` occupied (it must be empty) and update the neighbor counts around it.
  void occupy(Site s) noexcept;
  /// Mark `s` empty (it must be occupied) and update the neighbor counts around it.
  void vacate(Site s) noexcept;

 private:
  static constexpr std::uint8_t kOccupiedBit = 0x80;
  static constexpr std::uint8_t kCountMask = 0x0F;

  std::size_t occupied_count_{0};
  /// Per site: bit 7 = occupied, bits 0-3 = number of occupied neighbors.
  std::vector<std::uint8_t> state_;
};

}  // namespace rw
