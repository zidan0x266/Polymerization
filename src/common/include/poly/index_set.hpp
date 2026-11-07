// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// A set of indices with constant-time insertion, removal, membership test and access by
// position, so that a uniformly random member can be drawn in O(1).
#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace poly {

template <std::integral Id>
class IndexSet {
 public:
  IndexSet() = default;
  /// Empty set over the universe [0, universe).
  explicit IndexSet(std::size_t universe) : position_(universe, kAbsent) {}

  void reset(std::size_t universe) {
    members_.clear();
    position_.assign(universe, kAbsent);
  }

  [[nodiscard]] std::size_t size() const noexcept { return members_.size(); }
  [[nodiscard]] bool empty() const noexcept { return members_.empty(); }
  [[nodiscard]] bool contains(Id id) const noexcept { return position_[static_cast<std::size_t>(id)] != kAbsent; }
  /// Member at position `k` in [0, size()); positions change when members are removed.
  [[nodiscard]] Id operator[](std::size_t k) const noexcept { return members_[k]; }
  [[nodiscard]] const std::vector<Id>& members() const noexcept { return members_; }

  void insert(Id id) {
    if (contains(id)) return;
    position_[static_cast<std::size_t>(id)] = static_cast<std::uint32_t>(members_.size());
    members_.push_back(id);
  }
  void erase(Id id) noexcept {
    const std::uint32_t k = position_[static_cast<std::size_t>(id)];
    if (k == kAbsent) return;
    const Id last = members_.back();
    members_[k] = last;
    position_[static_cast<std::size_t>(last)] = k;
    members_.pop_back();
    position_[static_cast<std::size_t>(id)] = kAbsent;
  }
  void set(Id id, bool member) {
    if (member) {
      insert(id);
    } else {
      erase(id);
    }
  }

 private:
  static constexpr std::uint32_t kAbsent = std::numeric_limits<std::uint32_t>::max();
  std::vector<Id> members_;
  std::vector<std::uint32_t> position_;  ///< position in members_, or kAbsent
};

}  // namespace poly
