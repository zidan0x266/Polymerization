// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The reaction stage of a DLL time step, run after the motion stage. Reactions follow the shared
// FRP/ATRP chemistry (poly/chem/scheme.hpp): a unit bonds to a neighboring unit with the
// probability given by the rate constants, and a dormant initiator decomposes with k_init.
// Only pairs that include a reactive unit (a radical, or a dormant initiator in FRP) can react,
// so only these are attempted; the units that are reactive when the stage begins take part, in
// random order, and units activated during the stage wait for the next step.
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

#include "dll/model/system.hpp"
#include "poly/chem/progress.hpp"
#include "poly/chem/scheme.hpp"
#include "poly/index_set.hpp"
#include "poly/random.hpp"

namespace dll {

/// Which neighbors a unit attempts to react with in one time step.
enum class ReactionPartner {
  /// Every unit is the center of one attempt per step, with one of its twelve neighbors drawn at
  /// random: the rate constants are probabilities per attempt, as in the 3dbfm module.
  RandomNeighbor,
  /// Every pair of neighboring units is attempted once per step: the rate constants are
  /// probabilities per neighboring pair and time step.
  AllNeighbors
};

[[nodiscard]] std::string_view to_string(ReactionPartner p) noexcept;
[[nodiscard]] std::optional<ReactionPartner> parse_reaction_partner(std::string_view s) noexcept;

class Reactor {
 public:
  Reactor(const DllSystem& system, const poly::chem::ReactionScheme& scheme, ReactionPartner partner);

  /// One reaction stage. Returns the number of reactions.
  std::size_t step(DllSystem& system, poly::Rng& rng, poly::chem::ReactionProgress& progress);

  /// Number of units that can drive a reaction; zero means nothing can react any more.
  [[nodiscard]] std::size_t reactive_count() const noexcept { return reactive_.size(); }

 private:
  /// Attempt with `center` as the unit sampling and `partner` as the neighbor sampled.
  bool attempt_pair(DllSystem& system, Bead center, Bead partner, poly::Rng& rng, poly::chem::ReactionProgress& progress);
  bool attempt_initiation(DllSystem& system, Bead unit, poly::Rng& rng, poly::chem::ReactionProgress& progress);
  void step_random_neighbor(DllSystem& system, poly::Rng& rng, poly::chem::ReactionProgress& progress);
  void step_all_neighbors(DllSystem& system, poly::Rng& rng, poly::chem::ReactionProgress& progress);
  void update(const DllSystem& system, Bead unit);
  [[nodiscard]] bool in_stage(const DllSystem& system, Bead b) const noexcept;
  [[nodiscard]] Bead neighbor(const DllSystem& system, Bead b, int k) const noexcept;

  const poly::chem::ReactionScheme* scheme_;
  ReactionPartner partner_;
  poly::IndexSet<std::uint32_t> reactive_;  ///< by unit index (bead - first_unit)
  std::vector<std::uint32_t> order_;        ///< units taking part in the current stage
  std::vector<std::uint32_t> stage_of_;     ///< per unit: last stage it took part in
  std::uint32_t stage_{0};
  std::size_t reactions_{0};
};

}  // namespace dll
