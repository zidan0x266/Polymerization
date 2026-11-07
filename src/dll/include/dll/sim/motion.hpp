// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// One time step of the dynamic lattice liquid (Pakula; Polanowski & Pakula, J. Chem. Phys. 2002,
// 117, 4022). Every bead attempts a move to a random one of its twelve neighbors. With every site
// occupied, a move succeeds only as part of a closed loop of at least three beads, each taking the
// site of the next. Two beads trying to swap form no loop. Loops whose moves would leave a chain
// bond longer than a nearest-neighbor step are removed, which can break bonds in other loops, so
// the removal is repeated until no remaining loop breaks a bond. All remaining loops then move at
// once.
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "dll/model/system.hpp"
#include "poly/random.hpp"

namespace dll {

struct StepStats {
  std::size_t loops{0};     ///< loops that moved
  std::size_t rejected{0};  ///< loops removed because they would break a bond
  std::size_t moved{0};     ///< beads that moved
};

class Motion {
 public:
  explicit Motion(const DllSystem& system);

  /// One time step with a random field of attempted directions.
  StepStats step(DllSystem& system, poly::Rng& rng);
  /// One time step with the given direction (0..11, an index into kNeighbors) for every site.
  StepStats step(DllSystem& system, std::span<const std::uint8_t> directions);

 private:
  void find_targets(const DllSystem& system);
  void find_loops();
  void remove_bond_breaking_loops(const DllSystem& system, StepStats& stats);
  [[nodiscard]] bool breaks_bond(const DllSystem& system, std::uint32_t loop) const noexcept;
  [[nodiscard]] Vec3i new_position(const DllSystem& system, Bead b) const noexcept;
  [[nodiscard]] bool moves(Site s) const noexcept;

  std::vector<std::uint8_t> dir_;        ///< attempted direction per site
  std::vector<Site> target_;             ///< site each bead attempts to enter
  std::vector<std::uint32_t> mark_;      ///< loop search: walk that first visited a site
  std::vector<std::uint32_t> loop_of_;   ///< loop through a site, or kNoLoop
  std::vector<std::uint32_t> loop_begin_;
  std::vector<Site> loop_sites_;         ///< loops one after another, in moving order
  std::vector<std::uint8_t> rejected_;   ///< per loop
  std::vector<std::uint32_t> check_;
  std::vector<std::uint32_t> violators_;
  std::vector<std::uint32_t> queued_;    ///< per loop: round in which it was last queued
};

}  // namespace dll
