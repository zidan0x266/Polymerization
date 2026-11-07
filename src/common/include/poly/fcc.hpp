// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Geometry of a periodic face-centered cubic (FCC) lattice, shared by the lattice modules.
// Sites are stored in doubled integer coordinates: a site is (x, y, z) with 0 <= x, y, z < 2n and
// x + y + z even, where n is the number of cubic cells per edge. The twelve nearest neighbors
// are the permutations of (±1, ±1, 0). With nearest-neighbor distance b, the cubic cell edge is
// a = √2 b and a site sits at (b/√2)(x, y, z).
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <stdexcept>

namespace poly {

/// Small integer vector for lattice positions and steps (doubled coordinates).
struct Vec3i {
  int x{0};
  int y{0};
  int z{0};

  constexpr bool operator==(const Vec3i&) const noexcept = default;
  constexpr Vec3i operator+(Vec3i o) const noexcept { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3i operator-(Vec3i o) const noexcept { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3i operator-() const noexcept { return {-x, -y, -z}; }
  constexpr Vec3i& operator+=(Vec3i o) noexcept {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  [[nodiscard]] constexpr int dot(Vec3i o) const noexcept { return x * o.x + y * o.y + z * o.z; }
  [[nodiscard]] constexpr int norm2() const noexcept { return dot(*this); }
};

/// Index of a lattice site, 0 <= site < 4n^3.
using Site = std::uint32_t;

class FccGeometry {
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

  /// A lattice of `cells` cubic cells per edge (kMinCells to kMaxCells).
  explicit FccGeometry(int cells) : cells_(cells), extent_(2 * cells) {
    if (cells < kMinCells || cells > kMaxCells) {
      throw std::invalid_argument(
          std::format("FCC lattice needs {} to {} cells per edge, got {}", kMinCells, kMaxCells, cells));
    }
  }

  [[nodiscard]] int cells() const noexcept { return cells_; }
  /// Period of the doubled coordinates (2n).
  [[nodiscard]] int extent() const noexcept { return extent_; }
  [[nodiscard]] std::size_t num_sites() const noexcept {
    return 4 * static_cast<std::size_t>(cells_) * static_cast<std::size_t>(cells_) * static_cast<std::size_t>(cells_);
  }

  /// Map a lattice point at most one period outside the box into it.
  [[nodiscard]] Vec3i wrap(Vec3i p) const noexcept { return {wrap1(p.x), wrap1(p.y), wrap1(p.z)}; }
  /// Map any lattice point into the box.
  [[nodiscard]] Vec3i wrap_any(Vec3i p) const noexcept { return {modulo(p.x), modulo(p.y), modulo(p.z)}; }
  /// Shortest periodic image of a difference vector, components in [-n, n).
  [[nodiscard]] Vec3i minimum_image(Vec3i d) const noexcept {
    const auto shortest = [this](int c) { return modulo(c + cells_) - cells_; };
    return {shortest(d.x), shortest(d.y), shortest(d.z)};
  }
  /// Site index of a wrapped lattice point.
  [[nodiscard]] Site site(Vec3i p) const noexcept {
    return static_cast<Site>((static_cast<std::size_t>(p.x) * extent_ + p.y) * cells_ + (p.z >> 1));
  }
  /// Wrapped lattice point of a site.
  [[nodiscard]] Vec3i position(Site s) const noexcept {
    const auto zz = static_cast<int>(s % static_cast<Site>(cells_));
    const Site column = s / static_cast<Site>(cells_);
    const auto y = static_cast<int>(column % static_cast<Site>(extent_));
    const auto x = static_cast<int>(column / static_cast<Site>(extent_));
    return {x, y, 2 * zz + ((x + y) & 1)};
  }
  /// Site next to the wrapped point `p` along neighbor vector `k`.
  [[nodiscard]] Site neighbor(Vec3i p, int k) const noexcept { return site(wrap(p + kNeighbors[k])); }

 private:
  [[nodiscard]] int wrap1(int c) const noexcept {
    if (c < 0) return c + extent_;
    if (c >= extent_) return c - extent_;
    return c;
  }
  [[nodiscard]] int modulo(int c) const noexcept { return ((c % extent_) + extent_) % extent_; }

  int cells_;
  int extent_;
};

}  // namespace poly
