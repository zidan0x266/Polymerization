// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
#pragma once

#include <cstdint>

#include "bfm/core/types.hpp"
#include "bfm/model/unit.hpp"

namespace bfm {

enum class EventKind : std::uint8_t {
  Propagation,  ///< a radical adds across a vinyl group: a new bond, one vinyl consumed
  Termination,  ///< two radicals recombine: a new bond, no vinyl consumed
  Initiation    ///< a dormant initiator becomes a primary radical (unimolecular)
};

/// A reaction chosen by a ReactionScheme, to be executed by the simulation on the system.
struct ReactionEvent {
  EventKind kind{EventKind::Propagation};
  /// Propagation: the unit donating its radical. Termination: either radical. Initiation: the initiator.
  Index first{kNone};
  /// Propagation: the unit whose vinyl is attacked. Termination: the other radical. Initiation: kNone.
  Index second{kNone};
  /// One-based reaction category for the per-category event histogram (legacy `rcateg` layout).
  int category{0};
};

/// What a ReactionScheme is allowed to know about a unit when deciding.
struct UnitView {
  UnitState state;
  std::uint8_t nbonds{0};
  [[nodiscard]] constexpr int legacy_code() const noexcept { return legacy_type_code(state, nbonds); }
};

}  // namespace bfm
