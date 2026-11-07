// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// A set of unit indices with constant-time insertion, removal, membership test and access by
// position, so that a uniformly random member can be drawn in O(1).
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "bfm/core/types.hpp"

namespace bfm {

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
  [[nodiscard]] bool contains(Index id) const noexcept { return position_[static_cast<std::size_t>(id)] != kAbsent; }
  /// Member at position `k` in [0, size()); positions change when members are removed.
  [[nodiscard]] Index operator[](std::size_t k) const noexcept { return members_[k]; }

  void insert(Index id) {
    if (contains(id)) return;
    position_[static_cast<std::size_t>(id)] = static_cast<std::int32_t>(members_.size());
    members_.push_back(id);
  }
  void erase(Index id) noexcept {
    const std::int32_t k = position_[static_cast<std::size_t>(id)];
    if (k == kAbsent) return;
    const Index last = members_.back();
    members_[static_cast<std::size_t>(k)] = last;
    position_[static_cast<std::size_t>(last)] = k;
    members_.pop_back();
    position_[static_cast<std::size_t>(id)] = kAbsent;
  }
  void set(Index id, bool member) {
    if (member) {
      insert(id);
    } else {
      erase(id);
    }
  }

 private:
  static constexpr std::int32_t kAbsent = -1;
  std::vector<Index> members_;
  std::vector<std::int32_t> position_;  ///< position in members_, or kAbsent
};

}  // namespace bfm
