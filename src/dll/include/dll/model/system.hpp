// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The dynamic lattice liquid (DLL) system: a periodic FCC lattice on which every site holds
// exactly one bead. A bead is a unit of an inert precursor chain, a reactive unit (initiator,
// monomer or cross-linker, with its chemical state) or solvent. Bonds join beads on
// nearest-neighbor sites; reactions add bonds between units. Positions are also kept
// unwrapped, so displacements and periodic images are exact.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "poly/chem/progress.hpp"
#include "poly/chem/tables.hpp"
#include "poly/chem/unit.hpp"
#include "poly/fcc.hpp"

namespace dll {

using poly::FccGeometry;
using poly::Site;
using poly::Vec3i;
using poly::chem::ConnectivityTally;
using poly::chem::StructuralCounts;
using poly::chem::UnitState;
using poly::chem::UnitView;

/// Bead index. Beads are numbered chain after chain, then the reactive units, then the solvent
/// in site order.
using Bead = std::uint32_t;
inline constexpr int kMaxBonds = 4;

enum class BeadType : std::uint8_t { Solvent = 0, Chain = 1, Unit = 2 };

/// A bond between two beads, by zero-based bead index.
using BondPair = std::pair<Bead, Bead>;

class DllSystem {
 public:
  /// Chains in solvent: `chains` holds the unwrapped lattice points of the chain beads,
  /// `chain_length` per chain, joined by backbone bonds; every other site gets a solvent bead.
  DllSystem(int cells, std::size_t chain_length, std::span<const Vec3i> chains);
  /// Chains as above, plus reactive units in the states `units` at the unwrapped lattice points
  /// `unit_points`, plus the `bonds` formed by reactions (each listed once). Throws if two beads
  /// share a site, a bond does not join neighbors, or a bead would get more than kMaxBonds bonds.
  DllSystem(int cells, std::size_t chain_length, std::span<const Vec3i> chains, std::span<const Vec3i> unit_points,
            std::span<const UnitState> units, std::span<const BondPair> bonds = {});

  [[nodiscard]] const FccGeometry& lattice() const noexcept { return lattice_; }
  /// Number of beads, equal to the number of sites.
  [[nodiscard]] std::size_t size() const noexcept { return site_.size(); }
  [[nodiscard]] std::size_t chain_length() const noexcept { return chain_length_; }
  [[nodiscard]] std::size_t n_chains() const noexcept { return n_chains_; }
  [[nodiscard]] std::size_t chain_beads() const noexcept { return chain_length_ * n_chains_; }
  [[nodiscard]] std::size_t n_units() const noexcept { return states_.size(); }
  /// Chain beads and units: every bead that is not solvent.
  [[nodiscard]] std::size_t solute_beads() const noexcept { return chain_beads() + n_units(); }
  /// Bead index of the first reactive unit.
  [[nodiscard]] Bead first_unit() const noexcept { return static_cast<Bead>(chain_beads()); }

  [[nodiscard]] Bead occupant(Site s) const noexcept { return occupant_[s]; }
  [[nodiscard]] Site site(Bead b) const noexcept { return site_[b]; }
  [[nodiscard]] Vec3i unwrapped(Bead b) const noexcept { return unwrapped_[b]; }
  [[nodiscard]] BeadType type(Bead b) const noexcept {
    if (b < chain_beads()) return BeadType::Chain;
    return b < solute_beads() ? BeadType::Unit : BeadType::Solvent;
  }
  [[nodiscard]] std::span<const Bead> bonds(Bead b) const noexcept { return {bonds_[b].data(), nbonds_[b]}; }
  [[nodiscard]] bool bonded(Bead a, Bead b) const noexcept;

  /// Chemical state of the unit `b` (type(b) must be Unit).
  [[nodiscard]] UnitState& state(Bead b) noexcept { return states_[b - first_unit()]; }
  [[nodiscard]] const UnitState& state(Bead b) const noexcept { return states_[b - first_unit()]; }
  [[nodiscard]] UnitView view(Bead b) const noexcept { return {state(b), nbonds_[b]}; }
  [[nodiscard]] std::span<const UnitState> states() const noexcept { return states_; }

  /// Bond the units `a` and `b`, which must sit on neighboring sites and have a free bond slot,
  /// and count the new bond in the connectivity tally.
  void bond(Bead a, Bead b);
  /// A representative bead of the molecule that `b` (not solvent) belongs to.
  [[nodiscard]] Bead molecule(Bead b) const noexcept;
  [[nodiscard]] const ConnectivityTally& connectivity() const noexcept { return tally_; }
  void reset_connectivity_tally() noexcept { tally_ = {}; }
  /// Number of units carrying each legacy structural-unit code.
  [[nodiscard]] StructuralCounts structural_counts() const noexcept;

  /// Shift the beads along a closed loop of sites: the bead on loop[k] moves to loop[k + 1]
  /// (cyclically) by the nearest-neighbor step kNeighbors[directions[loop[k]]].
  void shift_loop(std::span<const Site> loop, std::span<const std::uint8_t> directions) noexcept;

  /// Check that every site holds one bead, positions match sites, bonds join neighbors and
  /// the unit states are possible.
  [[nodiscard]] std::expected<void, std::string> validate() const;

 private:
  void place(Bead b, Vec3i p);
  void link(Bead a, Bead b);
  void unite(Bead a, Bead b) noexcept;

  FccGeometry lattice_;
  std::size_t chain_length_;
  std::size_t n_chains_;
  std::vector<Bead> occupant_;  ///< per site
  std::vector<Site> site_;      ///< per bead
  std::vector<Vec3i> unwrapped_;
  std::vector<std::array<Bead, kMaxBonds>> bonds_;
  std::vector<std::uint8_t> nbonds_;
  std::vector<UnitState> states_;  ///< per unit
  std::vector<Bead> parent_;       ///< molecules: union-find forest over chain beads and units
  std::vector<std::uint32_t> molecule_size_;
  ConnectivityTally tally_{};
};

}  // namespace dll
