// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// One time step of the dynamic lattice liquid (Pakula; Polanowski & Pakula, J. Chem. Phys. 2002,
// 117, 4022). Every bead attempts a move to a random one of its twelve neighbors. With every site
// occupied, a move succeeds only as part of a closed loop of at least three beads, each taking the
// site of the next. Two beads trying to swap form no loop. A loop is removed if its moves would
// leave a bond longer than a nearest-neighbor step, or would make chains pass through each other:
// two bonds touching, or a bead without bonds passing through a bond, while all beads move in
// straight lines at equal speed (see crossing.hpp). Loops that break a bond are removed first, in
// rounds, since removing a loop changes the moves around it; then loops that would let chains
// cross, given the loops left; and so on until no remaining loop fails. All remaining loops then
// move at once.
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "dll/model/system.hpp"
#include "dll/sim/crossing.hpp"
#include "poly/random.hpp"

namespace dll {

struct StepStats {
  std::size_t loops{0};     ///< loops that moved
  std::size_t rejected{0};  ///< loops removed because they would break a bond
  std::size_t blocked{0};   ///< loops removed because chains would pass through each other
  std::size_t moved{0};     ///< beads that moved
};

class Motion {
 public:
  /// `allow_crossing` lets chains pass through each other (phantom chains), for comparison only.
  explicit Motion(const DllSystem& system, bool allow_crossing = false);

  /// One time step with a random field of attempted directions.
  StepStats step(DllSystem& system, poly::Rng& rng);
  /// One time step with the given direction (0..11, an index into kNeighbors) for every site.
  StepStats step(DllSystem& system, std::span<const std::uint8_t> directions);

 private:
  /// A bead with its position relative to the bond being checked, which unfolds the periodic box.
  struct Located {
    Bead bead;
    Vec3i start;
  };

  void find_targets(const DllSystem& system);
  void find_loops();
  void remove_failed_loops(const DllSystem& system, StepStats& stats);
  void remove_bond_breaking_loops(const DllSystem& system);
  void find_crossings(const DllSystem& system, bool all);
  void reject_flagged(const DllSystem& system);
  [[nodiscard]] bool breaks_bond(const DllSystem& system, std::uint32_t loop) const noexcept;
  void check_moving_bond(const DllSystem& system, Bead a, Bead b);
  void check_changed_bond(const DllSystem& system, Bead y, Bead z);
  template <class Visit>
  void for_each_nearby(const DllSystem& system, Located x, Located y, Visit&& visit) const;
  [[nodiscard]] Vec3i new_position(const DllSystem& system, Bead b) const noexcept;
  [[nodiscard]] bool moves(Site s) const noexcept;
  [[nodiscard]] Vec3i step_of(const DllSystem& system, Bead b) const noexcept;
  [[nodiscard]] Path path(const DllSystem& system, Located x) const noexcept { return {x.start, step_of(system, x.bead)}; }
  /// Bead on site `s` after the step, given the loops accepted so far.
  [[nodiscard]] Bead arrival(const DllSystem& system, Site s) const noexcept;
  void flag(std::uint32_t loop, bool crossing);
  void flag_bead(const DllSystem& system, Bead b);
  [[nodiscard]] bool tested_before(Bead y, Bead z, Vec3i y_start);

  bool allow_crossing_;
  std::vector<std::uint8_t> dir_;        ///< attempted direction per site
  std::vector<Site> target_;             ///< site each bead attempts to enter
  std::vector<Site> pred_;               ///< for a site on a loop: the site whose bead enters it
  std::vector<std::uint32_t> mark_;      ///< loop search: walk that first visited a site
  std::vector<std::uint32_t> loop_of_;   ///< loop through a site, or kNoLoop
  std::vector<std::uint32_t> loop_begin_;
  std::vector<Site> loop_sites_;         ///< loops one after another, in moving order
  std::vector<std::uint8_t> rejected_;   ///< per loop
  std::vector<std::uint32_t> check_;
  std::vector<std::uint32_t> violators_;
  std::vector<std::uint32_t> queued_;    ///< per loop: round in which it was last queued
  std::vector<std::uint32_t> flagged_;   ///< per loop: round in which it was last flagged
  std::vector<Bead> changed_;            ///< beads with bonds whose loop was removed since the last crossing check
  std::vector<std::uint8_t> involved_;   ///< per site: its bead moves or is bonded to one that does
  std::vector<Site> involved_sites_;
  struct Tested {
    Bead y, z;
    Vec3i y_start;
  };
  std::vector<Tested> tested_;  ///< bonds already tested against the current one
  std::uint32_t round_{0};
  std::size_t bond_failures_{0}, crossing_failures_{0};
};

}  // namespace dll
