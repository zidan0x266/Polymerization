// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/model/fcc.hpp"

#include <cassert>

namespace rw {

void FccLattice::occupy(Site s) noexcept {
  assert(!occupied(s));
  state_[s] |= kOccupiedBit;
  ++occupied_count_;
  const Vec3i p = position(s);
  for (int k = 0; k < kCoordination; ++k) ++state_[neighbor(p, k)];
}

void FccLattice::vacate(Site s) noexcept {
  assert(occupied(s));
  state_[s] &= static_cast<std::uint8_t>(~kOccupiedBit);
  --occupied_count_;
  const Vec3i p = position(s);
  for (int k = 0; k < kCoordination; ++k) --state_[neighbor(p, k)];
}

}  // namespace rw
