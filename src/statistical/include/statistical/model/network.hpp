// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The network of the statistical model: units with their chemical states and bonds, without
// coordinates. Molecules are tracked with the union-find of src/common; their representative
// units also count the pendant vinyl groups and the radicals of the molecule.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <vector>

#include "poly/chem/progress.hpp"
#include "poly/chem/tables.hpp"
#include "poly/chem/unit.hpp"
#include "poly/topology/bond_table.hpp"
#include "poly/topology/molecules.hpp"

namespace statistical {

using poly::chem::ConnectivityTally;
using poly::chem::StructuralCounts;
using poly::chem::UnitState;
using poly::chem::UnitView;

/// Unit index, 0 <= unit < size().
using Unit = std::uint32_t;
inline constexpr int kMaxBonds = 4;

/// Pendant vinyl groups of a unit: vinyls left on a unit that has already reacted.
[[nodiscard]] constexpr int pendant_vinyls(UnitState s) noexcept { return s.has_pendant_vinyl() ? s.vinyls : 0; }

class Network {
 public:
  explicit Network(std::vector<UnitState> units);

  [[nodiscard]] std::size_t size() const noexcept { return states_.size(); }
  [[nodiscard]] const UnitState& state(Unit u) const noexcept { return states_[u]; }
  [[nodiscard]] UnitView view(Unit u) const noexcept { return {states_[u], bonds_.count(u)}; }
  [[nodiscard]] std::span<const UnitState> states() const noexcept { return states_; }
  [[nodiscard]] std::span<const Unit> bonds(Unit u) const noexcept { return bonds_.partners(u); }
  [[nodiscard]] bool bonded(Unit a, Unit b) const noexcept { return bonds_.bonded(a, b); }

  /// Change the chemical state of `u`, keeping the counts of its molecule.
  void set_state(Unit u, UnitState s) noexcept;
  /// Bond `a` and `b` (distinct, not yet bonded, each with a free bond slot) and count the new
  /// bond in the connectivity tally.
  void bond(Unit a, Unit b);

  /// Representative unit of the molecule of `u`.
  [[nodiscard]] Unit molecule(Unit u) const noexcept { return molecules_.find(u); }
  /// Pendant vinyl groups and radicals on the molecule of `u`.
  [[nodiscard]] std::uint32_t molecule_pendants(Unit u) const noexcept { return pendants_[molecule(u)]; }
  [[nodiscard]] std::uint32_t molecule_radicals(Unit u) const noexcept { return radicals_[molecule(u)]; }

  [[nodiscard]] const ConnectivityTally& connectivity() const noexcept { return tally_; }
  void reset_connectivity_tally() noexcept { tally_ = {}; }
  /// Number of units carrying each legacy structural-unit code.
  [[nodiscard]] StructuralCounts structural_counts() const noexcept;

  /// Check bonds, molecules and the counts kept at the roots.
  [[nodiscard]] std::expected<void, std::string> validate() const;

 private:
  std::vector<UnitState> states_;
  poly::topology::BondTable<Unit, kMaxBonds> bonds_;
  poly::topology::Molecules<Unit> molecules_;
  std::vector<std::uint32_t> pendants_;  ///< at representatives: pendant vinyl groups in the molecule
  std::vector<std::uint32_t> radicals_;  ///< at representatives: radicals in the molecule
  ConnectivityTally tally_{};
};

}  // namespace statistical
