// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The dynamic lattice liquid (DLL) system: a periodic FCC lattice on which every site holds
// exactly one bead. Beads are chain units joined by bonds between nearest-neighbor sites, or
// solvent. Positions are also kept unwrapped, so displacements and periodic images are exact.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <vector>

#include "poly/fcc.hpp"

namespace dll {

using poly::FccGeometry;
using poly::Site;
using poly::Vec3i;

/// Bead index; beads are numbered chain after chain, then the solvent in site order.
using Bead = std::uint32_t;
inline constexpr int kMaxBonds = 4;

enum class BeadType : std::uint8_t { Solvent = 0, Chain = 1 };

class DllSystem {
 public:
  /// `chains` holds the unwrapped lattice points of the chain beads, `chain_length` per chain;
  /// every other site gets a solvent bead. Throws if two beads share a site or a bond is not a
  /// nearest-neighbor step.
  DllSystem(int cells, std::size_t chain_length, std::span<const Vec3i> chains);

  [[nodiscard]] const FccGeometry& lattice() const noexcept { return lattice_; }
  /// Number of beads, equal to the number of sites.
  [[nodiscard]] std::size_t size() const noexcept { return site_.size(); }
  [[nodiscard]] std::size_t chain_length() const noexcept { return chain_length_; }
  [[nodiscard]] std::size_t n_chains() const noexcept { return n_chains_; }
  [[nodiscard]] std::size_t chain_beads() const noexcept { return chain_length_ * n_chains_; }

  [[nodiscard]] Bead occupant(Site s) const noexcept { return occupant_[s]; }
  [[nodiscard]] Site site(Bead b) const noexcept { return site_[b]; }
  [[nodiscard]] Vec3i unwrapped(Bead b) const noexcept { return unwrapped_[b]; }
  [[nodiscard]] BeadType type(Bead b) const noexcept { return b < chain_beads() ? BeadType::Chain : BeadType::Solvent; }
  [[nodiscard]] std::span<const Bead> bonds(Bead b) const noexcept { return {bonds_[b].data(), nbonds_[b]}; }

  /// Shift the beads along a closed loop of sites: the bead on loop[k] moves to loop[k + 1]
  /// (cyclically) by the nearest-neighbor step kNeighbors[directions[loop[k]]].
  void shift_loop(std::span<const Site> loop, std::span<const std::uint8_t> directions) noexcept;

  /// Check that every site holds one bead, positions match sites and bonds join neighbors.
  [[nodiscard]] std::expected<void, std::string> validate() const;

 private:
  void add_bond(Bead a, Bead b) noexcept;

  FccGeometry lattice_;
  std::size_t chain_length_;
  std::size_t n_chains_;
  std::vector<Bead> occupant_;  ///< per site
  std::vector<Site> site_;      ///< per bead
  std::vector<Vec3i> unwrapped_;
  std::vector<std::array<Bead, kMaxBonds>> bonds_;
  std::vector<std::uint8_t> nbonds_;
};

}  // namespace dll
