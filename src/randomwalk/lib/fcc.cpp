// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/lattice/fcc.hpp"

#include <cassert>
#include <format>
#include <stdexcept>

namespace rw {

FccLattice::FccLattice(int cells) : cells_(cells), extent_(2 * cells) {
  if (cells < kMinCells || cells > kMaxCells) {
    throw std::invalid_argument(std::format("FCC lattice needs {} to {} cells per edge, got {}", kMinCells, kMaxCells, cells));
  }
  state_.assign(4 * static_cast<std::size_t>(cells) * cells * cells, 0);
}

Vec3i FccLattice::position(Site s) const noexcept {
  const auto zz = static_cast<int>(s % static_cast<Site>(cells_));
  const Site column = s / static_cast<Site>(cells_);
  const auto y = static_cast<int>(column % static_cast<Site>(extent_));
  const auto x = static_cast<int>(column / static_cast<Site>(extent_));
  return {x, y, 2 * zz + ((x + y) & 1)};
}

void FccLattice::occupy(Site s) noexcept {
  assert(!occupied(s));
  state_[s] |= kOccupiedBit;
  ++occupied_count_;
  const Vec3i p = position(s);
  for (const Vec3i& d : kNeighbors) ++state_[site(step(p, d))];
}

void FccLattice::vacate(Site s) noexcept {
  assert(occupied(s));
  state_[s] &= static_cast<std::uint8_t>(~kOccupiedBit);
  --occupied_count_;
  const Vec3i p = position(s);
  for (const Vec3i& d : kNeighbors) --state_[site(step(p, d))];
}

}  // namespace rw
