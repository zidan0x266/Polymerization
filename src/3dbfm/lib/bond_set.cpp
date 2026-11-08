// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/model/bond_set.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace bfm {

BondVectorSet::BondVectorSet() {
  generate_legacy_order();
  pair_antiparallel_vectors();
  build_lookup();
  build_shift_table();
  build_reaction_offsets();
}

// Reproduces legacy `bdibfl()`: for each of the six base vectors, emit 2 x 3 x 8 variants
// (two index swaps x three cyclic rotations x eight sign patterns) in a fixed order, keep the
// first occurrence of each distinct vector, append to the global list.
void BondVectorSet::generate_legacy_order() {
  constexpr std::array<Vec3i, 6> base{{{2, 0, 0}, {2, 1, 0}, {2, 1, 1}, {2, 2, 1}, {3, 0, 0}, {3, 1, 0}}};
  int filled = 0;
  for (Vec3i start : base) {
    std::vector<Vec3i> variants;
    variants.reserve(48);
    for (int swap = 0; swap < 2; ++swap) {
      for (int rotation = 0; rotation < 3; ++rotation) {
        const int a = start.x, b = start.y, c = start.z;
        variants.push_back({a, b, c});
        variants.push_back({a, b, -c});
        variants.push_back({a, -b, c});
        variants.push_back({a, -b, -c});
        variants.push_back({-a, b, c});
        variants.push_back({-a, b, -c});
        variants.push_back({-a, -b, c});
        variants.push_back({-a, -b, -c});
        start = {start.y, start.z, start.x};  // cyclic rotation
      }
      start = {start.y, start.x, start.z};  // transposition of the first two components
    }
    std::vector<Vec3i> distinct;
    for (const Vec3i& v : variants) {
      if (std::find(distinct.begin(), distinct.end(), v) == distinct.end()) distinct.push_back(v);
    }
    for (const Vec3i& v : distinct) {
      if (filled >= kNumBondVectors) throw std::logic_error("bond vector generation overflow");
      vectors_[filled++] = v;
    }
  }
  if (filled != kNumBondVectors) throw std::logic_error("bond vector generation produced != 108 vectors");
}

// Reproduces legacy `iniffxyz()`: arrange the set so that vector(107 - i) == -vector(i).
void BondVectorSet::pair_antiparallel_vectors() {
  for (int i = 0; i < kNumBondVectors / 2; ++i) {
    for (int j = 0; j < kNumBondVectors; ++j) {
      if (vectors_[j] == -vectors_[i]) std::swap(vectors_[j], vectors_[kNumBondVectors - 1 - i]);
    }
  }
  for (int i = 0; i < kNumBondVectors; ++i) {
    if (vectors_[kNumBondVectors - 1 - i] != -vectors_[i]) {
      throw std::logic_error("antiparallel pairing of bond vectors failed");
    }
  }
}

void BondVectorSet::build_lookup() {
  lookup_.fill(-1);
  for (int i = 0; i < kNumBondVectors; ++i) lookup_[lookup_slot(vectors_[i])] = static_cast<std::int16_t>(i);
}

std::optional<BondIndex> BondVectorSet::index_of(Vec3i v) const noexcept {
  if (!in_lookup_range(v)) return std::nullopt;
  const auto slot = lookup_[lookup_slot(v)];
  if (slot < 0) return std::nullopt;
  return static_cast<BondIndex>(slot);
}

void BondVectorSet::build_shift_table() {
  for (int i = 0; i < kNumBondVectors; ++i) {
    for (int d = 0; d < kNumDirections; ++d) {
      const auto target = index_of(vectors_[i] + unit_step(static_cast<Direction>(d)));
      shift_table_[i][d] = target ? static_cast<std::int16_t>(*target) : std::int16_t{-1};
    }
  }
}

// Reproduces legacy `ffxyz()`: lexicographic enumeration of (dx,dy,dz) in [-2,2]^3 with
// 4 <= |v|^2 <= 6. The stored bond index is that of -offset, i.e. the vector pointing from the
// partner anchored at `center + offset` back to the center unit.
void BondVectorSet::build_reaction_offsets() {
  int n = 0;
  for (int dx = -2; dx <= 2; ++dx) {
    for (int dy = -2; dy <= 2; ++dy) {
      for (int dz = -2; dz <= 2; ++dz) {
        const Vec3i offset{dx, dy, dz};
        const int l2 = offset.norm2();
        if (l2 < 4 || l2 > 6) continue;
        const auto bond = index_of(-offset);
        if (!bond) throw std::logic_error("reaction offset is not an allowed bond vector");
        if (n >= kNumReactionOffsets) throw std::logic_error("too many reaction offsets");
        reaction_offsets_[n++] = ReactionOffset{offset, *bond};
      }
    }
  }
  if (n != kNumReactionOffsets) throw std::logic_error("expected 54 reaction offsets");
  for (int i = 0; i < n; ++i) {
    if (reaction_offsets_[opposite_offset(i)].offset != -reaction_offsets_[i].offset) {
      throw std::logic_error("reaction offsets are not ordered symmetrically");
    }
  }
}

}  // namespace bfm
