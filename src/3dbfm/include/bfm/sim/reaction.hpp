// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The reactions: the FRP/ATRP chemistry shared with the dll and statistical modules
// (poly/chem/scheme.hpp and poly/chem/progress.hpp), and the sampling of reaction attempts on
// the lattice (map_reaction_draw).
#pragma once

#include <cstdint>
#include <optional>

#include "bfm/model/system.hpp"
#include "poly/chem/progress.hpp"
#include "poly/chem/scheme.hpp"

namespace bfm {

using poly::chem::CompetitionRule;
using poly::chem::EventKind;
using poly::chem::make_scheme;
using poly::chem::parse_competition_rule;
using poly::chem::RadicalChainScheme;
using poly::chem::RateConstants;
using poly::chem::Reaction;
using poly::chem::ReactionProgress;
using poly::chem::ReactionScheme;
using poly::chem::to_string;


/// A legacy reaction attempt: the center unit and the reaction-offset slot taken from it.
struct ReactionAttempt {
  Index center{kNone};
  int slot{0};
  constexpr bool operator==(const ReactionAttempt&) const noexcept = default;
};

/// Map a uniform draw (x, slot), x in [0, n) and slot in [0, 54), to the legacy attempt
/// (center, slot) it stands for, or to nothing when no reactive unit can be involved. With A the
/// set of reactive units:
///   x <  |A|          center A[x] at `slot` (the legacy attempt from a reactive center);
///   |A| <= x < 2|A|   the reactive unit A[x - |A|] is the partner: the center is the unit at
///                     offset `slot` from it, taken only if it is not reactive itself (otherwise
///                     the forward case covers the pair);
///   x >= 2|A|         nothing.
/// Every legacy attempt that involves a reactive unit is produced by exactly one draw, so its
/// probability stays 1/(54 n). When 2|A| > n the draw is used as in the legacy code.
[[nodiscard]] std::optional<ReactionAttempt> map_reaction_draw(const BfmSystem& system, const IndexSet& reactive,
                                                               std::uint32_t x, int slot) noexcept;

}  // namespace bfm
