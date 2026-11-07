// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "statistical/model/network.hpp"

#include <algorithm>
#include <format>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace statistical {

Network::Network(std::vector<UnitState> units)
    : states_(std::move(units)),
      bonds_(states_.size()),
      nbonds_(states_.size(), 0),
      parent_(states_.size()),
      size_(states_.size(), 1),
      pendants_(states_.size()),
      radicals_(states_.size()) {
  std::iota(parent_.begin(), parent_.end(), Unit{0});
  for (std::size_t u = 0; u < states_.size(); ++u) {
    pendants_[u] = static_cast<std::uint32_t>(pendant_vinyls(states_[u]));
    radicals_[u] = states_[u].radicals;
  }
}

bool Network::bonded(Unit a, Unit b) const noexcept { return std::ranges::find(bonds(a), b) != bonds(a).end(); }

Unit Network::molecule(Unit u) const noexcept {
  while (parent_[u] != u) u = parent_[u];
  return u;
}

void Network::set_state(Unit u, UnitState s) noexcept {
  const Unit root = molecule(u);
  pendants_[root] = pendants_[root] - static_cast<std::uint32_t>(pendant_vinyls(states_[u])) + static_cast<std::uint32_t>(pendant_vinyls(s));
  radicals_[root] = radicals_[root] - states_[u].radicals + s.radicals;
  states_[u] = s;
}

void Network::bond(Unit a, Unit b) {
  if (a == b || bonded(a, b) || nbonds_[a] >= kMaxBonds || nbonds_[b] >= kMaxBonds) {
    throw std::logic_error(std::format("cannot bond units {} and {}", a + 1, b + 1));
  }
  Unit ra = molecule(a), rb = molecule(b);
  if (nbonds_[a] == 0 && nbonds_[b] == 0) {
    ++tally_.monomer_monomer;
  } else if (nbonds_[a] != 0 && nbonds_[b] != 0) {
    ++(ra == rb ? tally_.intramolecular : tally_.polymer_polymer);
  } else {
    ++tally_.polymer_monomer;
  }
  bonds_[a][nbonds_[a]++] = b;
  bonds_[b][nbonds_[b]++] = a;
  if (ra == rb) return;
  // Union by size keeps the trees shallow, so molecule() needs no path compression.
  if (size_[ra] < size_[rb]) std::swap(ra, rb);
  parent_[rb] = ra;
  size_[ra] += size_[rb];
  pendants_[ra] += pendants_[rb];
  radicals_[ra] += radicals_[rb];
}

StructuralCounts Network::structural_counts() const noexcept {
  StructuralCounts counts{};
  for (std::size_t u = 0; u < states_.size(); ++u) {
    const int code = poly::chem::legacy_type_code(states_[u], nbonds_[u]);
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
    if (parent_[u] != u) continue;
    if (units[u] != size_[u] || pendants[u] != pendants_[u] || radicals[u] != radicals_[u]) {
      return fail(std::format("the counts of molecule {} are wrong", u + 1));
    }
  }
  return {};
}

}  // namespace statistical
