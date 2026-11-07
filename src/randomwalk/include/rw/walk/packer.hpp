// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Packing of precursor chains and cross-linkers on the FCC lattice. Each chain is a
// self-avoiding random walk with restricted bond angles started from a random empty site; the
// cross-linkers then take random empty sites. The original program is the Amorphous Polymer
// Generator of M. A. Tschopp and D. K. Ward, Mississippi State University (Hossain et al.,
// Polymer 2010, 51, 6071).
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "poly/random.hpp"
#include "rw/lattice/fcc.hpp"
#include "rw/model/packing.hpp"

namespace rw {

/// How the next site of a walk is chosen among the allowed empty neighbors of its head.
enum class SiteWeighting {
  /// As in the original program: weight 12 minus the number of neighboring beads from which a
  /// chain has already grown. The weights are taken, in the program's scan order, from the
  /// list of all empty neighbors but applied to the shorter list of those with an allowed
  /// angle, so a weight moves to another candidate whenever the angle rule removes a site.
  Legacy,
  /// With probability proportional to the number of empty neighbors of the candidate site,
  /// which steers chains away from crowded regions and dead ends.
  FreeNeighbors,
  /// With equal probability.
  Uniform
};

[[nodiscard]] std::string_view to_string(SiteWeighting w) noexcept;
[[nodiscard]] std::optional<SiteWeighting> parse_site_weighting(std::string_view s) noexcept;

struct WalkRules {
  /// The original program meant to allow 90 and 120 degrees, but its floating-point test
  /// "angle > 60" computes a 60-degree fold as 60.00000007 and lets it through.
  AngleSet angles{60, 90, 120};
  SiteWeighting weighting{SiteWeighting::Legacy};
  /// Attempts per chain; a chain whose head is trapped is removed and grown again elsewhere.
  int max_attempts{1000};
};

/// Order in which the original program scanned the twelve neighbors of a head at lattice point
/// `head`, as indices into FccLattice::kNeighbors. Exposed for testing.
[[nodiscard]] const std::array<int, FccLattice::kCoordination>& legacy_scan_order(Vec3i head) noexcept;

/// Smallest number of cells per edge for which `chain_beads` fill less than `occupancy` of the
/// lattice sites, i.e. floor(cbrt(N / (4 occupancy))) + 1 as in the original program.
[[nodiscard]] int cells_for_occupancy(std::size_t chain_beads, double occupancy);

class Packer {
 public:
  Packer(int cells, WalkRules rules, poly::Rng& rng);

  /// Grow `n_chains` chains of `chain_length` beads, then place `n_crosslinkers` single
  /// beads. Throws std::runtime_error if a chain cannot be grown in `max_attempts` attempts.
  [[nodiscard]] Packing pack(std::size_t n_chains, std::size_t chain_length, std::size_t n_crosslinkers,
                             double bond_length);

  /// One attempt to grow a chain of `beads.size()` beads; on a dead end the partial chain is
  /// removed again and false is returned.
  [[nodiscard]] bool grow_chain(std::span<Vec3i> beads);
  /// A uniformly random empty site.
  [[nodiscard]] Site random_free_site();

  [[nodiscard]] const FccLattice& lattice() const noexcept { return lattice_; }
  /// Number of chain attempts abandoned at a dead end so far.
  [[nodiscard]] std::size_t dead_ends() const noexcept { return dead_ends_; }

 private:
  /// Add `delta` to the head count of the twelve neighbors of `s` (legacy weighting).
  void count_head(Site s, int delta) noexcept;

  FccLattice lattice_;
  WalkRules rules_;
  poly::Rng& rng_;
  std::vector<Site> path_;  ///< sites of the chain being grown
  /// Legacy weighting only: per site, the number of neighboring beads that served as a head.
  std::vector<std::uint8_t> head_count_;
  std::size_t dead_ends_{0};
};

}  // namespace rw
