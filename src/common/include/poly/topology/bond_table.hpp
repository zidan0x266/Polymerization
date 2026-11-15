// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The bonds of every unit, at most MaxBonds each, kept inline so that walking the partners of a
// unit touches one small array. Every bond is listed from both ends.
#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace poly::topology {

template <std::integral Id, int MaxBonds>
class BondTable {
  static_assert(MaxBonds > 0 && MaxBonds < 256, "bond counts are stored in one byte");

 public:
  static constexpr int kMaxBonds = MaxBonds;

  BondTable() = default;
  /// `n` units without bonds.
  explicit BondTable(std::size_t n) : partners_(n), count_(n, 0) {}

  [[nodiscard]] std::size_t size() const noexcept { return count_.size(); }
  [[nodiscard]] std::span<const Id> partners(Id u) const noexcept { return {partners_[at(u)].data(), count_[at(u)]}; }
  [[nodiscard]] std::uint8_t count(Id u) const noexcept { return count_[at(u)]; }
  [[nodiscard]] bool full(Id u) const noexcept { return count_[at(u)] >= MaxBonds; }
  [[nodiscard]] bool bonded(Id a, Id b) const noexcept { return std::ranges::find(partners(a), b) != partners(a).end(); }

  /// Add the bond a-b to both lists; both units must have a free slot.
  void link(Id a, Id b) noexcept {
    partners_[at(a)][count_[at(a)]++] = b;
    partners_[at(b)][count_[at(b)]++] = a;
  }

 private:
  [[nodiscard]] static constexpr std::size_t at(Id u) noexcept { return static_cast<std::size_t>(u); }

  std::vector<std::array<Id, MaxBonds>> partners_;
  std::vector<std::uint8_t> count_;
};

}  // namespace poly::topology
