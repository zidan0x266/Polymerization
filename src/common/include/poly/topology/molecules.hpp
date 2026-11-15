// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Molecules as disjoint sets of units (union-find), for telling inter- from intramolecular
// reactions. Union by size keeps the trees shallow, so finding the molecule of a unit needs no
// path compression and stays const. Modules that keep totals per molecule (pendant groups,
// radicals) store them at the representative units and merge them as unite() reports.
#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

namespace poly::topology {

template <std::integral Id>
class Molecules {
 public:
  /// The representatives of two molecules just joined: `kept` represents the joined molecule.
  struct Merge {
    Id kept;
    Id absorbed;
  };

  Molecules() = default;
  explicit Molecules(std::size_t n) { reset(n); }

  /// `n` units, each a molecule of its own.
  void reset(std::size_t n) {
    parent_.resize(n);
    std::iota(parent_.begin(), parent_.end(), Id{0});
    units_.assign(n, 1);
  }

  [[nodiscard]] std::size_t size() const noexcept { return parent_.size(); }
  /// Representative unit of the molecule of `u`.
  [[nodiscard]] Id find(Id u) const noexcept {
    while (parent_[at(u)] != u) u = parent_[at(u)];
    return u;
  }
  [[nodiscard]] bool same(Id a, Id b) const noexcept { return find(a) == find(b); }
  [[nodiscard]] bool is_representative(Id u) const noexcept { return parent_[at(u)] == u; }
  /// Number of units in the molecule represented by `root`.
  [[nodiscard]] std::uint32_t units(Id root) const noexcept { return units_[at(root)]; }

  /// Join the molecules of `a` and `b`; std::nullopt if they are one molecule already. The
  /// larger molecule keeps its representative, the molecule of `a` on a tie.
  std::optional<Merge> unite(Id a, Id b) noexcept {
    Id ra = find(a), rb = find(b);
    if (ra == rb) return std::nullopt;
    if (units_[at(ra)] < units_[at(rb)]) std::swap(ra, rb);
    parent_[at(rb)] = ra;
    units_[at(ra)] += units_[at(rb)];
    return Merge{ra, rb};
  }

 private:
  [[nodiscard]] static constexpr std::size_t at(Id u) noexcept { return static_cast<std::size_t>(u); }

  std::vector<Id> parent_;
  std::vector<std::uint32_t> units_;  ///< at representatives: units in the molecule
};

}  // namespace poly::topology
