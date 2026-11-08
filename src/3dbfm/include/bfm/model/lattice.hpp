// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Periodic occupancy lattice for the bond fluctuation model. Every unit occupies the 2x2x2
// cube whose lowest corner ("anchor") is stored by the owning Monomer. Two arrays describe the
// lattice:
//   - one occupancy bit per site, for the excluded-volume tests of every move attempt; at one
//     bit per site a 100^3 lattice takes 125 KB and stays in the fastest cache levels;
//   - the index of the unit anchored at each site (kNone elsewhere), for finding reaction
//     partners and checking consistency.
#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bfm/model/types.hpp"

namespace bfm {

class PeriodicLattice {
 public:
  /// Largest supported extent: anchors are stored as 16-bit integers (see Monomer).
  static constexpr int kMaxExtent = 32766;

  PeriodicLattice(int lx, int ly, int lz) : dim_{lx, ly, lz} {
    for (int d : dim_) {
      if (d < 4 || d % 2 != 0 || d > kMaxExtent) {
        throw std::invalid_argument("lattice dimensions must be even, at least 4 and at most 32766");
      }
    }
    stride_y_ = static_cast<std::size_t>(lz);
    stride_x_ = static_cast<std::size_t>(ly) * stride_y_;
    const std::size_t sites = static_cast<std::size_t>(lx) * stride_x_;
    occupied_.assign((sites + 63) / 64, 0);
    anchors_.assign(sites, kNone);
  }

  [[nodiscard]] int lx() const noexcept { return dim_[0]; }
  [[nodiscard]] int ly() const noexcept { return dim_[1]; }
  [[nodiscard]] int lz() const noexcept { return dim_[2]; }
  [[nodiscard]] int extent(int axis) const noexcept { return dim_[axis]; }
  [[nodiscard]] std::size_t site_count() const noexcept { return anchors_.size(); }
  [[nodiscard]] Vec3i dimensions() const noexcept { return {dim_[0], dim_[1], dim_[2]}; }

  /// Linear index of a site (z runs fastest).
  [[nodiscard]] std::size_t site_index(int x, int y, int z) const noexcept {
    return static_cast<std::size_t>(x) * stride_x_ + static_cast<std::size_t>(y) * stride_y_ + static_cast<std::size_t>(z);
  }
  [[nodiscard]] std::size_t site_index(Vec3i p) const noexcept { return site_index(p.x, p.y, p.z); }

  [[nodiscard]] bool occupied(int x, int y, int z) const noexcept { return test(site_index(x, y, z)); }
  [[nodiscard]] bool occupied(Vec3i p) const noexcept { return test(site_index(p)); }
  /// C++23 multidimensional subscript: `lattice[x, y, z]` is true for an occupied site.
  [[nodiscard]] bool operator[](int x, int y, int z) const noexcept { return occupied(x, y, z); }

  /// The unit whose cube is anchored exactly at the (wrapped) site `p`, or kNone.
  [[nodiscard]] Index anchored_unit(Vec3i p) const noexcept { return anchors_[site_index(p)]; }

  /// The unit whose cube is anchored exactly at wrap(anchor + offset), or kNone.
  [[nodiscard]] Index unit_anchored_at(Vec3i anchor, Vec3i offset) const noexcept {
    return anchored_unit(wrap(anchor + offset));
  }

  /// Number of occupied sites (eight per unit in a consistent system).
  [[nodiscard]] std::size_t occupied_count() const noexcept {
    std::size_t n = 0;
    for (std::uint64_t w : occupied_) n += static_cast<std::size_t>(std::popcount(w));
    return n;
  }
  /// Number of sites at which some unit is anchored (one per unit in a consistent system).
  [[nodiscard]] std::size_t anchor_count() const noexcept {
    return static_cast<std::size_t>(std::ranges::count_if(anchors_, [](Index id) { return id != kNone; }));
  }

  /// Wrap a coordinate that lies within one period of the box on either side.
  [[nodiscard]] int wrap(int axis, int v) const noexcept {
    const int L = dim_[axis];
    if (v < 0) return v + L;
    if (v >= L) return v - L;
    return v;
  }
  [[nodiscard]] Vec3i wrap(Vec3i p) const noexcept { return {wrap(0, p.x), wrap(1, p.y), wrap(2, p.z)}; }
  [[nodiscard]] int next(int axis, int v) const noexcept { return v + 1 == dim_[axis] ? 0 : v + 1; }
  [[nodiscard]] int prev(int axis, int v) const noexcept { return v == 0 ? dim_[axis] - 1 : v - 1; }

  /// True if none of the eight sites of the cube anchored at `a` is occupied.
  [[nodiscard]] bool cube_is_free(Vec3i a) const noexcept {
    return all_free(x_face(a.x, a)) && all_free(x_face(next(0, a.x), a));
  }

  void occupy_cube(Vec3i a, Index id) noexcept {
    for (const Face& f : {x_face(a.x, a), x_face(next(0, a.x), a)}) {
      for (std::size_t s : f) set_bit(s);
    }
    anchors_[site_index(a)] = id;
  }
  void clear() noexcept {
    std::ranges::fill(occupied_, std::uint64_t{0});
    std::ranges::fill(anchors_, kNone);
  }

  /// Are the four sites the cube would newly cover when jumping in direction `d` all empty?
  [[nodiscard]] bool front_is_free(Vec3i a, Direction d) const noexcept { return all_free(entering_face(a, d)); }

  /// Move the cube of unit `id` anchored at `a` one site in direction `d` (the front sites
  /// must be free). Returns the new anchor.
  Vec3i shift_cube(Vec3i a, Direction d, Index id) noexcept {
    Face leaving{};
    Vec3i to = a;
    switch (d) {
      case Direction::PlusX: leaving = x_face(a.x, a); to.x = next(0, a.x); break;
      case Direction::MinusX: leaving = x_face(next(0, a.x), a); to.x = prev(0, a.x); break;
      case Direction::PlusY: leaving = y_face(a.y, a); to.y = next(1, a.y); break;
      case Direction::MinusY: leaving = y_face(next(1, a.y), a); to.y = prev(1, a.y); break;
      case Direction::PlusZ: leaving = z_face(a.z, a); to.z = next(2, a.z); break;
      case Direction::MinusZ: leaving = z_face(next(2, a.z), a); to.z = prev(2, a.z); break;
    }
    for (std::size_t s : entering_face(a, d)) set_bit(s);
    for (std::size_t s : leaving) clear_bit(s);
    anchors_[site_index(a)] = kNone;
    anchors_[site_index(to)] = id;
    return to;
  }

 private:
  /// Four sites of one face of a cube.
  using Face = std::array<std::size_t, 4>;

  /// The four sites of the cube anchored at `a` that lie in the plane x = `x` (resp. y, z).
  [[nodiscard]] Face x_face(int x, Vec3i a) const noexcept {
    const int y1 = next(1, a.y), z1 = next(2, a.z);
    return {site_index(x, a.y, a.z), site_index(x, y1, a.z), site_index(x, a.y, z1), site_index(x, y1, z1)};
  }
  [[nodiscard]] Face y_face(int y, Vec3i a) const noexcept {
    const int x1 = next(0, a.x), z1 = next(2, a.z);
    return {site_index(a.x, y, a.z), site_index(x1, y, a.z), site_index(a.x, y, z1), site_index(x1, y, z1)};
  }
  [[nodiscard]] Face z_face(int z, Vec3i a) const noexcept {
    const int x1 = next(0, a.x), y1 = next(1, a.y);
    return {site_index(a.x, a.y, z), site_index(x1, a.y, z), site_index(a.x, y1, z), site_index(x1, y1, z)};
  }
  /// The face a cube anchored at `a` would newly cover after a jump in direction `d`.
  [[nodiscard]] Face entering_face(Vec3i a, Direction d) const noexcept {
    switch (d) {
      case Direction::PlusX: return x_face(next(0, next(0, a.x)), a);
      case Direction::MinusX: return x_face(prev(0, a.x), a);
      case Direction::PlusY: return y_face(next(1, next(1, a.y)), a);
      case Direction::MinusY: return y_face(prev(1, a.y), a);
      case Direction::PlusZ: return z_face(next(2, next(2, a.z)), a);
      case Direction::MinusZ: return z_face(prev(2, a.z), a);
    }
    std::unreachable();
  }

  [[nodiscard]] bool test(std::size_t s) const noexcept { return ((occupied_[s >> 6] >> (s & 63)) & 1u) != 0; }
  /// Branch-free test of four sites: at high packing the outcome is hard to predict.
  [[nodiscard]] bool all_free(const Face& f) const noexcept {
    return (static_cast<unsigned>(test(f[0])) | static_cast<unsigned>(test(f[1])) | static_cast<unsigned>(test(f[2])) |
            static_cast<unsigned>(test(f[3]))) == 0;
  }
  void set_bit(std::size_t s) noexcept { occupied_[s >> 6] |= std::uint64_t{1} << (s & 63); }
  void clear_bit(std::size_t s) noexcept { occupied_[s >> 6] &= ~(std::uint64_t{1} << (s & 63)); }

  std::array<int, 3> dim_;
  std::size_t stride_x_{0};
  std::size_t stride_y_{0};
  std::vector<std::uint64_t> occupied_;  ///< one bit per site
  std::vector<Index> anchors_;           ///< unit anchored at each site, or kNone
};

}  // namespace bfm
