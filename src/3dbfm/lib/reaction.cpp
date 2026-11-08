// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/sim/reaction.hpp"

namespace bfm {

std::optional<ReactionAttempt> map_reaction_draw(const BfmSystem& system, const IndexSet& reactive, std::uint32_t x,
                                                 int slot) noexcept {
  const auto n = static_cast<std::uint32_t>(system.size());
  const auto na = static_cast<std::uint32_t>(reactive.size());
  if (2 * static_cast<std::uint64_t>(na) > n) return ReactionAttempt{static_cast<Index>(x), slot};
  if (x < na) return ReactionAttempt{reactive[x], slot};
  if (x < 2 * na) {
    const Index center = system.partner_at(reactive[x - na], slot);
    if (center == kNone || reactive.contains(center)) return std::nullopt;
    return ReactionAttempt{center, BondVectorSet::opposite_offset(slot)};
  }
  return std::nullopt;
}

}  // namespace bfm
