// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Fundamental value types shared by every module of the 3D bond-fluctuation-model
// (3DBFM) reactive Monte Carlo library.
#pragma once

#include <cstdint>
#include <string_view>
#include <utility>

namespace bfm {

/// Dense, zero-based monomer identifier. `kNone` marks an empty lattice site or "no partner".
using Index = std::int32_t;
inline constexpr Index kNone = -1;

/// Index into the ordered set of 108 allowed bond vectors (0..107).
using BondIndex = std::uint8_t;
inline constexpr int kNumBondVectors = 108;

/// Lattice offsets with |v|^2 in {4,5,6} at which a reaction partner may be anchored.
inline constexpr int kNumReactionOffsets = 54;

/// A unit takes part in at most four bonds (divinyl unit: two vinyl additions + two propagations).
inline constexpr int kMaxBondsPerUnit = 4;

/// Small integer vector used for lattice coordinates, offsets and bond vectors.
struct Vec3i {
  int x{0};
  int y{0};
  int z{0};

  constexpr bool operator==(const Vec3i&) const noexcept = default;
  constexpr Vec3i operator+(Vec3i o) const noexcept { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3i operator-(Vec3i o) const noexcept { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3i operator-() const noexcept { return {-x, -y, -z}; }
  constexpr Vec3i& operator+=(Vec3i o) noexcept {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  [[nodiscard]] constexpr int norm2() const noexcept { return x * x + y * y + z * z; }
  [[nodiscard]] constexpr int operator[](int axis) const noexcept {
    return axis == 0 ? x : (axis == 1 ? y : z);
  }
};

/// The six nearest-neighbor jump directions. The enumeration order equals the legacy
/// numbering (legacy 1..6 -> 0..5) so that tables and restart files stay comparable.
enum class Direction : std::uint8_t { PlusX = 0, PlusY = 1, PlusZ = 2, MinusZ = 3, MinusY = 4, MinusX = 5 };
inline constexpr int kNumDirections = 6;

[[nodiscard]] constexpr Vec3i unit_step(Direction d) noexcept {
  switch (d) {
    case Direction::PlusX: return {1, 0, 0};
    case Direction::PlusY: return {0, 1, 0};
    case Direction::PlusZ: return {0, 0, 1};
    case Direction::MinusZ: return {0, 0, -1};
    case Direction::MinusY: return {0, -1, 0};
    case Direction::MinusX: return {-1, 0, 0};
  }
  std::unreachable();
}

[[nodiscard]] constexpr std::string_view to_string(Direction d) noexcept {
  switch (d) {
    case Direction::PlusX: return "+x";
    case Direction::PlusY: return "+y";
    case Direction::PlusZ: return "+z";
    case Direction::MinusZ: return "-z";
    case Direction::MinusY: return "-y";
    case Direction::MinusX: return "-x";
  }
  std::unreachable();
}

}  // namespace bfm
