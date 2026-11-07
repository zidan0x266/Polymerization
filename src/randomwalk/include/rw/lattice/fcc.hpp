// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Periodic face-centered cubic (FCC) lattice. Sites are stored in doubled integer coordinates:
// a site is (x, y, z) with 0 <= x, y, z < 2n and x + y + z even, where n is the number of cubic
// cells per edge. The twelve nearest neighbors are the permutations of (±1, ±1, 0). With
// nearest-neighbor distance b (the bond length), the cubic cell edge is a = √2 b and a site sits
// at (b/√2)(x, y, z).
#pragma once

#include <array>
#include <cstdint>
#include <initializer_list>
#include <vector>

namespace rw {

/// Small integer vector for lattice positions and bond vectors (doubled coordinates).
struct Vec3i {
  int x{0};
  int y{0};
  int z{0};

  constexpr bool operator==(const Vec3i&) const noexcept = default;
  constexpr Vec3i operator+(Vec3i o) const noexcept { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3i operator-(Vec3i o) const noexcept { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3i operator-() const noexcept { return {-x, -y, -z}; }
  [[nodiscard]] constexpr int dot(Vec3i o) const noexcept { return x * o.x + y * o.y + z * o.z; }
  [[nodiscard]] constexpr int norm2() const noexcept { return dot(*this); }
};

/// Index of a lattice site, 0 <= site < 4n^3.
using Site = std::uint32_t;

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

class FccLattice {
 public:
  static constexpr int kCoordination = 12;
  /// Two cells keep the twelve neighbors of a site distinct; 4n^3 sites must fit a 32-bit index.
  static constexpr int kMinCells = 2;
  static constexpr int kMaxCells = 1000;
  /// Nearest-neighbor vectors, each with squared length 2.
  static constexpr std::array<Vec3i, kCoordination> kNeighbors{{{1, 1, 0},
                                                                {1, -1, 0},
                                                                {-1, 1, 0},
                                                                {-1, -1, 0},
                                                                {1, 0, 1},
                                                                {1, 0, -1},
                                                                {-1, 0, 1},
                                                                {-1, 0, -1},
                                                                {0, 1, 1},
                                                                {0, 1, -1},
                                                                {0, -1, 1},
                                                                {0, -1, -1}}};

  /// A lattice of `cells` cubic cells per edge (kMinCells to kMaxCells), all sites empty.
  explicit FccLattice(int cells);

  [[nodiscard]] int cells() const noexcept { return cells_; }
  /// Period of the doubled coordinates (2n).
  [[nodiscard]] int extent() const noexcept { return extent_; }
  [[nodiscard]] std::size_t num_sites() const noexcept { return state_.size(); }

  /// Map any lattice point (x + y + z even) into the periodic box.
  [[nodiscard]] Vec3i wrap(Vec3i p) const noexcept { return {wrap1(p.x), wrap1(p.y), wrap1(p.z)}; }
  /// Site index of a wrapped lattice point.
  [[nodiscard]] Site site(Vec3i p) const noexcept {
    return static_cast<Site>((static_cast<std::size_t>(p.x) * extent_ + p.y) * cells_ + (p.z >> 1));
  }
  /// Wrapped lattice point of a site.
  [[nodiscard]] Vec3i position(Site s) const noexcept;
  /// Site next to the wrapped point `p` along neighbor vector `k`.
  [[nodiscard]] Site neighbor(Vec3i p, int k) const noexcept { return site(step(p, kNeighbors[k])); }

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

  /// One wrapped coordinate, for points at most one period outside the box.
  [[nodiscard]] int wrap1(int c) const noexcept {
    if (c < 0) return c + extent_;
    if (c >= extent_) return c - extent_;
    return c;
  }
  /// Wrapped point `p + d` for a wrapped `p` and a neighbor vector `d`.
  [[nodiscard]] Vec3i step(Vec3i p, Vec3i d) const noexcept { return wrap(p + d); }

  int cells_;
  int extent_;
  std::size_t occupied_count_{0};
  /// Per site: bit 7 = occupied, bits 0-3 = number of occupied neighbors.
  std::vector<std::uint8_t> state_;
};

}  // namespace rw
