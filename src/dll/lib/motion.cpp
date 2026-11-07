// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/sim/motion.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <utility>

namespace dll {

namespace {

constexpr std::uint32_t kNoLoop = std::numeric_limits<std::uint32_t>::max();

/// Offsets from a site to itself, its twelve neighbors and its six second neighbors: every site
/// within √2 bond lengths.
constexpr std::array<Vec3i, 19> kBall = [] {
  std::array<Vec3i, 19> ball{};
  for (int k = 0; k < FccGeometry::kCoordination; ++k) ball[k + 1] = FccGeometry::kNeighbors[k];
  ball[13] = {2, 0, 0};
  ball[14] = {-2, 0, 0};
  ball[15] = {0, 2, 0};
  ball[16] = {0, -2, 0};
  ball[17] = {0, 0, 2};
  ball[18] = {0, 0, -2};
  return ball;
}();

}  // namespace

Motion::Motion(const DllSystem& system, bool allow_crossing)
    : allow_crossing_(allow_crossing),
      dir_(system.size()),
      target_(system.size()),
      pred_(system.size()),
      mark_(system.size()),
      loop_of_(system.size()),
      involved_(system.size(), 0) {}

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
  remove_failed_loops(system, stats);
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
    for (std::size_t i = begin; i < loop_sites_.size(); ++i) {
      loop_of_[loop_sites_[i]] = loop;
      pred_[target_[loop_sites_[i]]] = loop_sites_[i];
    }
    loop_begin_.push_back(static_cast<std::uint32_t>(loop_sites_.size()));
  }
  rejected_.assign(loop_begin_.size() - 1, 0);
}

bool Motion::moves(Site s) const noexcept {
  const std::uint32_t loop = loop_of_[s];
  return loop != kNoLoop && !rejected_[loop];
}

Vec3i Motion::step_of(const DllSystem& system, Bead b) const noexcept {
  const Site s = system.site(b);
  return moves(s) ? FccGeometry::kNeighbors[dir_[s]] : Vec3i{};
}

Vec3i Motion::new_position(const DllSystem& system, Bead b) const noexcept { return system.unwrapped(b) + step_of(system, b); }

Bead Motion::arrival(const DllSystem& system, Site s) const noexcept {
  return system.occupant(moves(s) ? pred_[s] : s);
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

void Motion::flag(std::uint32_t loop, bool crossing) {
  if (flagged_[loop] == round_) return;
  flagged_[loop] = round_;
  violators_.push_back(loop);
  ++(crossing ? crossing_failures_ : bond_failures_);
}

void Motion::flag_bead(const DllSystem& system, Bead b) {
  const Site s = system.site(b);
  if (moves(s)) flag(loop_of_[s], true);
}

bool Motion::tested_before(Bead y, Bead z, Vec3i y_start) {
  for (const Tested& t : tested_) {
    if (t.y == y && t.z == z && t.y_start == y_start) return true;
  }
  tested_.push_back({y, z, y_start});
  return false;
}

template <class Visit>
void Motion::for_each_nearby(const DllSystem& system, Located x, Located y, Visit&& visit) const {
  // Two segments of length at most one bond can touch only if an end of one comes within one
  // bond length of an end of the other at that instant. The ends move by at most one bond
  // length, so some pair of ends is then within √2 bond lengths at the start or at the end of
  // the step: the candidates are the beads on the sites around the ends before the step, and
  // those arriving on the sites around them after it. Sites whose bead neither moves nor is
  // bonded to a moving bead are skipped.
  const FccGeometry& lattice = system.lattice();
  for (const Located end : {x, y}) {
    const Vec3i before = lattice.position(system.site(end.bead));
    const Vec3i step = step_of(system, end.bead);
    const Vec3i after = lattice.wrap(before + step);
    for (const Vec3i o : kBall) {
      const Site s = lattice.site(lattice.wrap(before + o));
      if (involved_[s]) visit(Located{system.occupant(s), end.start + o});
    }
    for (const Vec3i o : kBall) {
      const Site s = lattice.site(lattice.wrap(after + o));
      if (!involved_[s]) continue;
      const Bead w = arrival(system, s);
      visit(Located{w, end.start + step + o - step_of(system, w)});
    }
  }
}

void Motion::check_moving_bond(const DllSystem& system, Bead a, Bead b) {
  // Only a bond whose two beads both move can touch another bond or a bead (shown by exhaustive
  // tests, tests/dll); the other one must move at least at one end.
  const FccGeometry& lattice = system.lattice();
  const Located la{a, {}};
  const Located lb{b, lattice.minimum_image(system.unwrapped(b) - system.unwrapped(a))};
  const Path pa = path(system, la), pb = path(system, lb);
  const std::pair key{std::min(a, b), std::max(a, b)};
  tested_.clear();
  for_each_nearby(system, la, lb, [&](Located w) {
    if (w.bead == a || w.bead == b) return;
    const bool w_moves = moves(system.site(w.bead));
    const auto partners = system.bonds(w.bead);
    if (partners.empty()) {  // a bead without bonds must not pass through the bond
      if (w_moves && !tested_before(w.bead, w.bead, w.start) && segments_touch(pa, pb, path(system, w), path(system, w))) {
        flag_bead(system, a);
        flag_bead(system, b);
        flag_bead(system, w.bead);
      }
      return;
    }
    for (const Bead v : partners) {
      if (v == a || v == b) continue;
      const bool v_moves = moves(system.site(v));
      if (!w_moves && !v_moves) continue;
      if (w_moves && v_moves && std::pair{std::min(w.bead, v), std::max(w.bead, v)} < key) continue;  // checked from that bond
      const Located lv{v, w.start + lattice.minimum_image(system.unwrapped(v) - system.unwrapped(w.bead))};
      const bool first = w.bead < v;
      if (tested_before(first ? w.bead : v, first ? v : w.bead, first ? w.start : lv.start)) continue;
      if (segments_touch(pa, pb, path(system, w), path(system, lv))) {
        for (const Bead q : {a, b, w.bead, v}) flag_bead(system, q);
      }
    }
  });
}

void Motion::check_changed_bond(const DllSystem& system, Bead y, Bead z) {
  // y no longer moves, z still does: the bond can now touch only bonds that move at both ends.
  const FccGeometry& lattice = system.lattice();
  const Located ly{y, {}};
  const Located lz{z, lattice.minimum_image(system.unwrapped(z) - system.unwrapped(y))};
  const Path py = path(system, ly), pz = path(system, lz);
  tested_.clear();
  for_each_nearby(system, ly, lz, [&](Located w) {
    if (w.bead == y || w.bead == z || !moves(system.site(w.bead))) return;
    for (const Bead v : system.bonds(w.bead)) {
      if (v == y || v == z || !moves(system.site(v))) continue;
      const Located lv{v, w.start + lattice.minimum_image(system.unwrapped(v) - system.unwrapped(w.bead))};
      const bool first = w.bead < v;
      if (tested_before(first ? w.bead : v, first ? v : w.bead, first ? w.start : lv.start)) continue;
      if (segments_touch(path(system, w), path(system, lv), py, pz)) {
        for (const Bead q : {z, w.bead, v}) flag_bead(system, q);
      }
    }
  });
}

void Motion::reject_flagged(const DllSystem& system) {
  // Remove the flagged loops; the loops bonded to them need another bond check, and the bonds
  // of their beads, which now stay, another crossing check.
  for (const std::uint32_t loop : violators_) rejected_[loop] = 1;
  check_.clear();
  for (const std::uint32_t loop : violators_) {
    for (std::uint32_t i = loop_begin_[loop]; i < loop_begin_[loop + 1]; ++i) {
      const Bead b = system.occupant(loop_sites_[i]);
      if (!system.bonds(b).empty()) changed_.push_back(b);
      for (const Bead q : system.bonds(b)) {
        const std::uint32_t other = loop_of_[system.site(q)];
        if (other == kNoLoop || rejected_[other] || queued_[other] == round_) continue;
        queued_[other] = round_;
        check_.push_back(other);
      }
    }
  }
  violators_.clear();
}

void Motion::remove_bond_breaking_loops(const DllSystem& system) {
  // Rounds of simultaneous removal: every loop that breaks a bond given the loops accepted at
  // the start of the round is removed; then only loops bonded to removed ones need another look.
  while (!check_.empty()) {
    ++round_;
    for (const std::uint32_t loop : check_) {
      if (!rejected_[loop] && breaks_bond(system, loop)) flag(loop, false);
    }
    reject_flagged(system);
  }
}

void Motion::find_crossings(const DllSystem& system, bool all) {
  ++round_;
  for (const Site s : involved_sites_) involved_[s] = 0;
  involved_sites_.clear();
  const auto involve = [&](Site s) {
    if (!involved_[s]) {
      involved_[s] = 1;
      involved_sites_.push_back(s);
    }
  };
  const auto n_loops = static_cast<std::uint32_t>(rejected_.size());
  for (std::uint32_t loop = 0; loop < n_loops; ++loop) {
    if (rejected_[loop]) continue;
    for (std::uint32_t i = loop_begin_[loop]; i < loop_begin_[loop + 1]; ++i) {
      involve(loop_sites_[i]);
      for (const Bead q : system.bonds(system.occupant(loop_sites_[i]))) involve(system.site(q));
    }
  }
  if (all) {
    for (std::uint32_t loop = 0; loop < n_loops; ++loop) {
      if (rejected_[loop]) continue;
      for (std::uint32_t i = loop_begin_[loop]; i < loop_begin_[loop + 1]; ++i) {
        const Bead a = system.occupant(loop_sites_[i]);
        for (const Bead b : system.bonds(a)) {
          if (a < b && moves(system.site(b))) check_moving_bond(system, a, b);
        }
      }
    }
  } else {
    for (const Bead y : changed_) {
      for (const Bead z : system.bonds(y)) {
        if (!moves(system.site(y)) && moves(system.site(z))) check_changed_bond(system, y, z);
      }
    }
  }
  changed_.clear();
}

void Motion::remove_failed_loops(const DllSystem& system, StepStats& stats) {
  const auto n_loops = static_cast<std::uint32_t>(rejected_.size());
  check_.resize(n_loops);
  std::iota(check_.begin(), check_.end(), 0u);
  queued_.assign(n_loops, 0);
  flagged_.assign(n_loops, 0);
  changed_.clear();
  violators_.clear();
  round_ = 0;
  bond_failures_ = crossing_failures_ = 0;
  remove_bond_breaking_loops(system);
  if (!allow_crossing_) {
    // The first crossing check covers every bond whose beads both move; after that only the
    // bonds changed by later removals can start to touch.
    for (bool all = true;; all = false) {
      find_crossings(system, all);
      if (violators_.empty()) break;
      reject_flagged(system);
      remove_bond_breaking_loops(system);
    }
  }
  stats.rejected += bond_failures_;
  stats.blocked += crossing_failures_;
}

}  // namespace dll
