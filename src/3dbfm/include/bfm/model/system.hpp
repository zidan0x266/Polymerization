// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The simulated system: all units, their connectivity, the occupancy lattice and the molecule
// bookkeeping needed to classify each new bond as monomer-monomer, polymer-monomer,
// polymer-polymer or intramolecular (the legacy chain_premark/mark12 tally). This class knows how
// to move units and how to bond them; it does not know any chemistry (see reaction/scheme.hpp).
#pragma once

#include <array>
#include <cstdint>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#include "bfm/core/random.hpp"
#include "bfm/core/types.hpp"
#include "bfm/lattice/bond_set.hpp"
#include "bfm/lattice/lattice.hpp"
#include "bfm/model/unit.hpp"

namespace bfm {

struct Bond {
  Index partner{kNone};
  /// Index of the bond vector pointing from `partner` to the owner of this record
  /// (the legacy `monbd` convention: moving the owner by +x shifts this vector by +x).
  BondIndex vector{0};
};

/// One unit, packed into 32 bytes so that two fit in a cache line: every move attempt touches
/// a random unit and, when it moves, its bonded partners. Periodic image counts, needed only
/// for output, are kept apart by BfmSystem.
struct alignas(32) Monomer {
  std::array<std::int16_t, 3> pos{};  ///< wrapped anchor of the 2x2x2 cube (zero-based)
  UnitState state;
  std::uint8_t nbonds{0};
  std::array<BondIndex, kMaxBondsPerUnit> vectors{};  ///< vectors[i]: from partners[i] to this unit
  std::array<Index, kMaxBondsPerUnit> partners{kNone, kNone, kNone, kNone};

  [[nodiscard]] constexpr Vec3i anchor() const noexcept { return {pos[0], pos[1], pos[2]}; }
  constexpr void set_anchor(Vec3i a) noexcept {
    pos = {static_cast<std::int16_t>(a.x), static_cast<std::int16_t>(a.y), static_cast<std::int16_t>(a.z)};
  }

  [[nodiscard]] constexpr Bond bond(int i) const noexcept { return {partners[i], vectors[i]}; }
  /// The bonds of this unit as a range of Bond values.
  [[nodiscard]] auto bonds() const noexcept {
    return std::views::iota(0, int{nbonds}) | std::views::transform([this](int i) { return bond(i); });
  }
  /// Position of `other` in the bond list, or -1.
  [[nodiscard]] constexpr int slot_of(Index other) const noexcept {
    for (int i = 0; i < nbonds; ++i) {
      if (partners[i] == other) return i;
    }
    return -1;
  }
  [[nodiscard]] constexpr bool bonded_to(Index other) const noexcept { return slot_of(other) >= 0; }
  constexpr void add_bond(Bond b) noexcept {
    partners[nbonds] = b.partner;
    vectors[nbonds] = b.vector;
    ++nbonds;
  }
};
static_assert(sizeof(Monomer) == 32, "Monomer is meant to fill exactly half a cache line");

/// Tally of how reactions connect pre-existing structures (legacy mon_mon/poly_mon/...).
struct ConnectivityTally {
  std::uint64_t monomer_monomer{0};
  std::uint64_t polymer_monomer{0};
  std::uint64_t polymer_polymer{0};
  std::uint64_t intramolecular{0};
};

class BfmSystem {
 public:
  BfmSystem(PeriodicLattice lattice, const BondVectorSet& bonds);

  // ---- construction -------------------------------------------------------------------------
  /// Place one cube per entry of `composition` at random, non-overlapping positions of the
  /// even sub-lattice (the legacy initial condition) and assign the states in order.
  void place_randomly(std::span<const UnitState> composition, Rng& rng);
  /// Rebuild from explicit monomers (restart). `unwrapped` holds the unwrapped anchors (empty:
  /// all units in the central image). Validates overlaps and bond consistency.
  void load(std::vector<Monomer> monomers, std::span<const Vec3i> unwrapped = {});

  // ---- Monte Carlo elementary moves ---------------------------------------------------------
  /// Try to jump unit `id` one lattice site in direction `d`. Returns true if accepted.
  bool attempt_move(Index id, Direction d) noexcept;

  /// Unit anchored exactly at the `slot`-th reaction offset from `id`, or kNone.
  [[nodiscard]] Index partner_at(Index id, int slot) const noexcept {
    return lattice_.unit_anchored_at(monomers_[id].anchor(), bonds_->reaction_offsets()[slot].offset);
  }

  /// Create a bond between `center` and the unit at reaction offset `slot` from it, updating
  /// connectivity tallies and the cluster structure. Both units must have a free bond slot.
  void bond_through_offset(Index center, Index partner, int slot);

  // ---- molecule bookkeeping for the reaction tally (union-find) ------------------------------
  [[nodiscard]] Index molecule_root(Index id) noexcept;
  [[nodiscard]] bool same_molecule(Index a, Index b) noexcept { return molecule_root(a) == molecule_root(b); }

  // ---- accessors ----------------------------------------------------------------------------
  [[nodiscard]] Index size() const noexcept { return static_cast<Index>(monomers_.size()); }
  [[nodiscard]] const Monomer& monomer(Index id) const noexcept { return monomers_[id]; }
  [[nodiscard]] Monomer& monomer(Index id) noexcept { return monomers_[id]; }
  [[nodiscard]] std::span<const Monomer> monomers() const noexcept { return monomers_; }
  /// Unwrapped anchor (wrapped anchor plus periodic images), for trajectories and restarts.
  [[nodiscard]] Vec3i unwrapped(Index id) const noexcept {
    const Vec3i a = monomers_[id].anchor();
    const Vec3i img = images_[id];
    return {a.x + img.x * lattice_.lx(), a.y + img.y * lattice_.ly(), a.z + img.z * lattice_.lz()};
  }
  [[nodiscard]] const PeriodicLattice& lattice() const noexcept { return lattice_; }
  [[nodiscard]] const BondVectorSet& bond_set() const noexcept { return *bonds_; }
  [[nodiscard]] const ConnectivityTally& connectivity() const noexcept { return tally_; }
  void reset_connectivity_tally() noexcept { tally_ = {}; }

  /// Number of units carrying each legacy structural-unit code 0..9 (the legacy s_0..s_9).
  [[nodiscard]] std::array<std::int64_t, kNumLegacyTypeCodes> structural_counts() const noexcept;

  /// Expensive consistency check (lattice occupancy and anchors, bond reciprocity and geometry).
  void validate() const;

 private:
  void add_bond(Index a, Index b, BondIndex vector_from_b_to_a);
  void unite(Index a, Index b) noexcept;

  PeriodicLattice lattice_;
  const BondVectorSet* bonds_;
  std::vector<Monomer> monomers_;
  std::vector<Vec3i> images_;  ///< periodic image of each unit; touched only when a move wraps
  std::vector<Index> parent_;
  ConnectivityTally tally_{};
};

}  // namespace bfm
