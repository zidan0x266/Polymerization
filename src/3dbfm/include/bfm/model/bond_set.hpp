// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The allowed bond-vector set of the 3D bond fluctuation model (Carmesin & Kremer; Deutsch &
// Binder): all 108 vectors obtained from {[2,0,0],[2,1,0],[2,1,1],[2,2,1],[3,0,0],[3,1,0]} by
// permutation and sign change. The generation order reproduces the legacy `bdibfl()` +
// `iniffxyz()` routines exactly, so that runs follow the legacy programs step by step.
#pragma once

#include <array>
#include <optional>
#include <span>
#include <utility>

#include "bfm/model/types.hpp"

namespace bfm {

class BondVectorSet {
 public:
  /// One of the 54 lattice offsets at which a reaction partner may be anchored, together with
  /// the index of the bond vector pointing from that partner back to the center unit.
  struct ReactionOffset {
    Vec3i offset;
    BondIndex bond_from_partner_to_center;
  };

  BondVectorSet();

  [[nodiscard]] const Vec3i& vector(BondIndex i) const noexcept { return vectors_[i]; }

  /// The antiparallel partner of a bond index: vector(reverse(i)) == -vector(i).
  [[nodiscard]] static constexpr BondIndex reverse(BondIndex i) noexcept {
    return static_cast<BondIndex>(kNumBondVectors - 1 - i);
  }

  /// Index of a vector in the set, if it belongs to it.
  [[nodiscard]] std::optional<BondIndex> index_of(Vec3i v) const noexcept;

  /// Index of the vector obtained by shifting bond `i` by one lattice unit in direction `d`,
  /// if the result is still an allowed bond vector (the legacy `move[][]` table).
  [[nodiscard]] std::optional<BondIndex> shifted(BondIndex i, Direction d) const noexcept {
    const auto v = shift_table_[i][std::to_underlying(d)];
    if (v < 0) return std::nullopt;
    return static_cast<BondIndex>(v);
  }

  [[nodiscard]] std::span<const ReactionOffset> reaction_offsets() const noexcept {
    return reaction_offsets_;
  }

  /// The slot of the opposite offset: reaction_offsets()[opposite_offset(s)].offset is
  /// -reaction_offsets()[s].offset (the offsets are generated in a sign-symmetric order).
  [[nodiscard]] static constexpr int opposite_offset(int slot) noexcept { return kNumReactionOffsets - 1 - slot; }

  /// All 108 vectors in canonical (legacy) order.
  [[nodiscard]] std::span<const Vec3i> vectors() const noexcept { return vectors_; }

 private:
  static constexpr int kLookupExtent = 7;  // coordinates in [-3, 3]
  [[nodiscard]] static constexpr int lookup_slot(Vec3i v) noexcept {
    return ((v.x + 3) * kLookupExtent + (v.y + 3)) * kLookupExtent + (v.z + 3);
  }
  [[nodiscard]] static constexpr bool in_lookup_range(Vec3i v) noexcept {
    return v.x >= -3 && v.x <= 3 && v.y >= -3 && v.y <= 3 && v.z >= -3 && v.z <= 3;
  }

  void generate_legacy_order();
  void pair_antiparallel_vectors();
  void build_lookup();
  void build_shift_table();
  void build_reaction_offsets();

  std::array<Vec3i, kNumBondVectors> vectors_{};
  std::array<std::int16_t, kLookupExtent * kLookupExtent * kLookupExtent> lookup_{};
  std::array<std::array<std::int16_t, kNumDirections>, kNumBondVectors> shift_table_{};
  std::array<ReactionOffset, kNumReactionOffsets> reaction_offsets_{};
};

}  // namespace bfm
