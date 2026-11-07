// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Geometry of chain crossing. During a time step every moving bead travels in a straight line
// from its site to the next one, all at the same uniform speed, and a bond is the straight
// segment between its two beads. Two bonds cross if the segments touch at some instant; a bead
// passes through a bond if it touches the segment. The test is exact up to rounding: contacts
// can only happen when the four ends are coplanar, so only the real roots of a cubic in time
// are examined (or, if the ends stay coplanar throughout, the instants at which an end meets
// the other segment).
#pragma once

#include "poly/fcc.hpp"

namespace dll {

using poly::Vec3i;

/// Straight path of a bead during a time step, in lattice coordinates.
struct Path {
  Vec3i start;
  Vec3i step;  ///< zero if the bead stays
  constexpr bool operator==(const Path&) const noexcept = default;
};

/// Whether the segment a–b touches the segment c–d at some instant of the step. Pass c == d
/// for a single bead. The four beads must be distinct and must not collide.
[[nodiscard]] bool segments_touch(Path a, Path b, Path c, Path d) noexcept;

}  // namespace dll
