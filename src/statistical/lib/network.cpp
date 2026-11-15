// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "statistical/model/network.hpp"

#include <format>
#include <stdexcept>
#include <utility>

namespace statistical {

Network::Network(std::vector<UnitState> units)
    : states_(std::move(units)),
      bonds_(states_.size()),
      molecules_(states_.size()),
      pendants_(states_.size()),
      radicals_(states_.size()) {
  for (std::size_t u = 0; u < states_.size(); ++u) {
    pendants_[u] = static_cast<std::uint32_t>(pendant_vinyls(states_[u]));
    radicals_[u] = states_[u].radicals;
  }
}

void Network::set_state(Unit u, UnitState s) noexcept {
  const Unit root = molecule(u);
  pendants_[root] = pendants_[root] - static_cast<std::uint32_t>(pendant_vinyls(states_[u])) + static_cast<std::uint32_t>(pendant_vinyls(s));
  radicals_[root] = radicals_[root] - states_[u].radicals + s.radicals;
  states_[u] = s;
}

void Network::bond(Unit a, Unit b) {
  if (a == b || bonded(a, b) || bonds_.full(a) || bonds_.full(b)) {
    throw std::logic_error(std::format("cannot bond units {} and {}", a + 1, b + 1));
  }
  tally_.count(bonds_.count(a), bonds_.count(b), [&] { return molecules_.same(a, b); });
  bonds_.link(a, b);
  if (const auto merge = molecules_.unite(a, b)) {
    pendants_[merge->kept] += pendants_[merge->absorbed];
    radicals_[merge->kept] += radicals_[merge->absorbed];
  }
}

StructuralCounts Network::structural_counts() const noexcept {
  StructuralCounts counts{};
  for (std::size_t u = 0; u < states_.size(); ++u) {
    const int code = poly::chem::legacy_type_code(states_[u], bonds_.count(static_cast<Unit>(u)));
    if (code >= 0 && code < poly::chem::kNumLegacyTypeCodes) ++counts[static_cast<std::size_t>(code)];
  }
  return counts;
}

std::expected<void, std::string> Network::validate() const {
  const auto fail = [](std::string message) { return std::unexpected(std::move(message)); };
  const std::size_t n = size();
  std::vector<std::uint32_t> units(n, 0), pendants(n, 0), radicals(n, 0);
  for (Unit u = 0; u < n; ++u) {
    const UnitState s = states_[u];
    if (s.vinyls > poly::chem::max_vinyls(s.species) || s.radicals > 2) return fail(std::format("unit {} has an impossible state", u + 1));
    for (const Unit p : bonds(u)) {
      if (p >= n || p == u || !bonded(p, u)) return fail(std::format("bond {}-{} is invalid or one-sided", u + 1, p + 1));
      if (molecule(p) != molecule(u)) return fail(std::format("bond {}-{} joins two molecules", u + 1, p + 1));
    }
    const Unit root = molecule(u);
    ++units[root];
    pendants[root] += static_cast<std::uint32_t>(pendant_vinyls(s));
    radicals[root] += s.radicals;
  }
  for (Unit u = 0; u < n; ++u) {
    if (!molecules_.is_representative(u)) continue;
    if (units[u] != molecules_.units(u) || pendants[u] != pendants_[u] || radicals[u] != radicals_[u]) {
      return fail(std::format("the counts of molecule {} are wrong", u + 1));
    }
  }
  return {};
}

}  // namespace statistical
