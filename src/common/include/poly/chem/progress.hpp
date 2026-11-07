// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Running totals of a polymerization and the state changes of an accepted reaction, the same in
// every lattice module. The modules form the bonds; this file keeps the chemical bookkeeping.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "poly/chem/scheme.hpp"
#include "poly/chem/unit.hpp"

namespace poly::chem {

/// How new bonds connect pre-existing structures (legacy mon_mon/poly_mon/...).
struct ConnectivityTally {
  std::uint64_t monomer_monomer{0};
  std::uint64_t polymer_monomer{0};
  std::uint64_t polymer_polymer{0};
  std::uint64_t intramolecular{0};
};

/// Running totals of the reaction phase for one sample.
struct ReactionProgress {
  std::uint64_t total_vinyls{0};     ///< N_mono + 2 N_di + N_initiator (legacy vinyl budget)
  std::uint64_t consumed_vinyls{0};  ///< counts the initiator-derived units as consumed (legacy definition)
  std::uint64_t propagations{0};
  std::uint64_t terminations{0};
  std::uint64_t initiations{0};
  std::uint64_t live_radicals{0};
  std::vector<std::uint64_t> category_counts;  ///< one-based categories; index 0 unused

  /// Start from zero with `n_categories` reaction categories.
  void reset(std::size_t n_categories) {
    *this = ReactionProgress{};
    category_counts.assign(n_categories + 1, 0);
  }

  /// Add a unit present at the start of the reaction phase. Initiator-derived units count as one
  /// consumed vinyl each (legacy budget), so the conversion starts at N_initiator / total.
  void count_unit(UnitState s) noexcept {
    const auto full = static_cast<std::uint64_t>(max_vinyls(s.species));
    total_vinyls += full;
    consumed_vinyls += full - s.vinyls;
    live_radicals += s.radicals;
    if (s.species == Species::Initiator) {
      ++total_vinyls;
      ++consumed_vinyls;
    }
  }

  [[nodiscard]] double conversion() const noexcept {
    return total_vinyls ? static_cast<double>(consumed_vinyls) / static_cast<double>(total_vinyls) : 0.0;
  }
};

/// Carry out the pair reaction `r` proposed for (a, b): update both states and the totals. The
/// caller forms the bond between the two units.
inline void apply(const Reaction& r, UnitState& a, UnitState& b, ReactionProgress& progress) noexcept {
  if (r.kind == EventKind::Propagation) {
    UnitState& attacker = r.second_attacks ? b : a;
    UnitState& target = r.second_attacks ? a : b;
    --attacker.radicals;
    --target.vinyls;
    ++target.radicals;
    ++progress.propagations;
    ++progress.consumed_vinyls;
  } else if (r.kind == EventKind::Termination) {
    --a.radicals;
    --b.radicals;
    ++progress.terminations;
    progress.live_radicals -= 2;
  }
  if (r.category > 0 && static_cast<std::size_t>(r.category) < progress.category_counts.size()) {
    ++progress.category_counts[static_cast<std::size_t>(r.category)];
  }
}

/// Carry out the unimolecular reaction `r` (initiation) of unit `a`.
inline void apply(const Reaction& r, UnitState& a, ReactionProgress& progress) noexcept {
  if (r.kind == EventKind::Initiation) {
    ++a.radicals;
    ++progress.initiations;
    ++progress.live_radicals;
  }
  if (r.category > 0 && static_cast<std::size_t>(r.category) < progress.category_counts.size()) {
    ++progress.category_counts[static_cast<std::size_t>(r.category)];
  }
}

}  // namespace poly::chem
