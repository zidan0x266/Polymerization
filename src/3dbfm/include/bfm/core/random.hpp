// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Random number source: the shared xoshiro256++ generator (src/common) plus the lattice draw
// used by the move step.
#pragma once

#include "bfm/core/types.hpp"
#include "poly/random.hpp"

namespace bfm {

using poly::splitmix64;
using poly::Xoshiro256pp;

class Rng : public poly::Rng {
 public:
  using poly::Rng::Rng;

  [[nodiscard]] Direction random_direction() noexcept { return static_cast<Direction>(below(kNumDirections)); }
};

}  // namespace bfm
