// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/sim/motion.hpp"

#include <algorithm>
#include <limits>

namespace dll {

namespace {

constexpr std::uint32_t kNoLoop = std::numeric_limits<std::uint32_t>::max();

}  // namespace

Motion::Motion(const DllSystem& system)
    : dir_(system.size()), target_(system.size()), mark_(system.size()), loop_of_(system.size()) {}

StepStats Motion::step(DllSystem& system, poly::Rng& rng) {
  const std::size_t n = dir_.size();
  std::size_t s = 0;
  for (; s + 1 < n; s += 2) {
    const auto [a, b] = rng.below2(FccGeometry::kCoordination, FccGeometry::kCoordination);
    dir_[s] = static_cast<std::uint8_t>(a);
    dir_[s + 1] = static_cast<std::uint8_t>(b);
  }
  if (s < n) dir_[s] = static_cast<std::uint8_t>(rng.below(FccGeometry::kCoordination));
  return step(system, dir_);
}

StepStats Motion::step(DllSystem& system, std::span<const std::uint8_t> directions) {
  if (directions.data() != dir_.data()) std::ranges::copy(directions, dir_.begin());
  StepStats stats;
  find_targets(system);
  find_loops();
  remove_bond_breaking_loops(system, stats);
  for (std::uint32_t loop = 0; loop + 1 < loop_begin_.size(); ++loop) {
    if (rejected_[loop]) continue;
    const std::span<const Site> sites(loop_sites_.data() + loop_begin_[loop], loop_begin_[loop + 1] - loop_begin_[loop]);
    system.shift_loop(sites, dir_);
    ++stats.loops;
    stats.moved += sites.size();
  }
  return stats;
}

void Motion::find_targets(const DllSystem& system) {
  // Sites are numbered x-major, then y, then z (see FccGeometry::site), so the positions can be
  // generated in order without divisions.
  const FccGeometry& lattice = system.lattice();
  const int extent = lattice.extent();
  const int cells = lattice.cells();
  Site s = 0;
  for (int x = 0; x < extent; ++x) {
    for (int y = 0; y < extent; ++y) {
      const int parity = (x + y) & 1;
      for (int zz = 0; zz < cells; ++zz, ++s) {
        target_[s] = lattice.neighbor({x, y, 2 * zz + parity}, dir_[s]);
      }
    }
  }
}

void Motion::find_loops() {
  // The attempts form a graph in which every site points to one neighbor; its cycles are the
  // closed loops. Following the pointers from every unvisited site finds each cycle once.
  std::ranges::fill(mark_, 0u);
  std::ranges::fill(loop_of_, kNoLoop);
  loop_begin_.assign(1, 0);
  loop_sites_.clear();
  const auto n = static_cast<Site>(target_.size());
  for (Site start = 0; start < n; ++start) {
    if (mark_[start] != 0) continue;
    const std::uint32_t walk = start + 1;
    Site u = start;
    while (mark_[u] == 0) {
      mark_[u] = walk;
      u = target_[u];
    }
    if (mark_[u] != walk) continue;  // ran into an earlier walk
    const std::size_t begin = loop_sites_.size();
    Site v = u;
    do {
      loop_sites_.push_back(v);
      v = target_[v];
    } while (v != u);
    if (loop_sites_.size() - begin < 3) {  // two beads trying to swap
      loop_sites_.resize(begin);
      continue;
    }
    const auto loop = static_cast<std::uint32_t>(loop_begin_.size() - 1);
    for (std::size_t i = begin; i < loop_sites_.size(); ++i) loop_of_[loop_sites_[i]] = loop;
    loop_begin_.push_back(static_cast<std::uint32_t>(loop_sites_.size()));
  }
  rejected_.assign(loop_begin_.size() - 1, 0);
}

bool Motion::moves(Site s) const noexcept {
  const std::uint32_t loop = loop_of_[s];
  return loop != kNoLoop && !rejected_[loop];
}

Vec3i Motion::new_position(const DllSystem& system, Bead b) const noexcept {
  const Site s = system.site(b);
  return moves(s) ? system.unwrapped(b) + FccGeometry::kNeighbors[dir_[s]] : system.unwrapped(b);
}

bool Motion::breaks_bond(const DllSystem& system, std::uint32_t loop) const noexcept {
  for (std::uint32_t i = loop_begin_[loop]; i < loop_begin_[loop + 1]; ++i) {
    const Bead b = system.occupant(loop_sites_[i]);
    const auto partners = system.bonds(b);
    if (partners.empty()) continue;
    const Vec3i p = new_position(system, b);
    for (const Bead q : partners) {
      if (system.lattice().minimum_image(new_position(system, q) - p).norm2() != 2) return true;
    }
  }
  return false;
}

void Motion::remove_bond_breaking_loops(const DllSystem& system, StepStats& stats) {
  // Rounds of simultaneous removal: every loop that breaks a bond given the loops accepted at the
  // start of the round is removed; then only loops bonded to removed ones need another look.
  const auto n_loops = static_cast<std::uint32_t>(rejected_.size());
  check_.resize(n_loops);
  for (std::uint32_t loop = 0; loop < n_loops; ++loop) check_[loop] = loop;
  queued_.assign(n_loops, 0);
  for (std::uint32_t round = 1; !check_.empty(); ++round) {
    violators_.clear();
    for (const std::uint32_t loop : check_) {
      if (!rejected_[loop] && breaks_bond(system, loop)) violators_.push_back(loop);
    }
    for (const std::uint32_t loop : violators_) rejected_[loop] = 1;
    stats.rejected += violators_.size();
    check_.clear();
    for (const std::uint32_t loop : violators_) {
      for (std::uint32_t i = loop_begin_[loop]; i < loop_begin_[loop + 1]; ++i) {
        for (const Bead q : system.bonds(system.occupant(loop_sites_[i]))) {
          const std::uint32_t other = loop_of_[system.site(q)];
          if (other == kNoLoop || rejected_[other] || queued_[other] == round) continue;
          queued_[other] = round;
          check_.push_back(other);
        }
      }
    }
  }
}

}  // namespace dll
