// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/model/system.hpp"

#include <format>
#include <limits>
#include <stdexcept>

namespace dll {

namespace {

constexpr Bead kNoBead = std::numeric_limits<Bead>::max();

}  // namespace

DllSystem::DllSystem(int cells, std::size_t chain_length, std::span<const Vec3i> chains)
    : DllSystem(cells, chain_length, chains, {}, {}, {}) {}

DllSystem::DllSystem(int cells, std::size_t chain_length, std::span<const Vec3i> chains, std::span<const Vec3i> unit_points,
                     std::span<const UnitState> units, std::span<const BondPair> bonds)
    : lattice_(cells),
      chain_length_(chain_length),
      n_chains_(chain_length ? chains.size() / chain_length : 0),
      states_(units.begin(), units.end()) {
  if (chain_length < 2 && !chains.empty()) throw std::invalid_argument("chains need at least two beads");
  if (chain_length && chains.size() % chain_length != 0) throw std::invalid_argument("bead count is not a multiple of the chain length");
  if (unit_points.size() != units.size()) throw std::invalid_argument("every unit needs one lattice point");
  const std::size_t n = lattice_.num_sites();
  if (chains.size() + units.size() > n) throw std::invalid_argument("more chain beads and units than lattice sites");

  occupant_.assign(n, kNoBead);
  site_.resize(n);
  unwrapped_.resize(n);
  bonds_ = poly::topology::BondTable<Bead, kMaxBonds>(n);
  molecules_.reset(solute_beads());

  for (Bead b = 0; b < chains.size(); ++b) {
    place(b, chains[b]);
    if (b % chain_length_ != 0) {
      if ((chains[b] - chains[b - 1]).norm2() != 2) {
        throw std::invalid_argument(std::format("bond {}-{} is not a nearest-neighbor step", b, b + 1));
      }
      link(b - 1, b);
      molecules_.unite(b - 1, b);
    }
  }
  for (std::size_t u = 0; u < units.size(); ++u) {
    const UnitState s = units[u];
    if (s.vinyls > poly::chem::max_vinyls(s.species) || s.radicals > 2) {
      throw std::invalid_argument(std::format("unit {} has an impossible chemical state", u + 1));
    }
    place(first_unit() + static_cast<Bead>(u), unit_points[u]);
  }
  for (const auto& [a, b] : bonds) {
    if (a >= solute_beads() || b >= solute_beads() || type(a) != BeadType::Unit || type(b) != BeadType::Unit || a == b) {
      throw std::invalid_argument(std::format("bond {}-{} does not join two units", a + 1, b + 1));
    }
    if (bonded(a, b)) throw std::invalid_argument(std::format("bond {}-{} is listed twice", a + 1, b + 1));
    if (lattice_.minimum_image(unwrapped_[b] - unwrapped_[a]).norm2() != 2) {
      throw std::invalid_argument(std::format("bond {}-{} does not join neighboring sites", a + 1, b + 1));
    }
    link(a, b);
    molecules_.unite(a, b);
  }

  auto solvent = static_cast<Bead>(solute_beads());
  for (Site s = 0; s < n; ++s) {
    if (occupant_[s] != kNoBead) continue;
    occupant_[s] = solvent;
    site_[solvent] = s;
    unwrapped_[solvent] = lattice_.position(s);
    ++solvent;
  }
}

void DllSystem::place(Bead b, Vec3i p) {
  const char* kind = b < chain_beads() ? "chain bead" : "unit";
  const Bead number = (b < chain_beads() ? b : b - first_unit()) + 1;
  if (((p.x + p.y + p.z) & 1) != 0) throw std::invalid_argument(std::format("{} {} is not on a lattice site", kind, number));
  const Site s = lattice_.site(lattice_.wrap_any(p));
  if (occupant_[s] != kNoBead) throw std::invalid_argument(std::format("{} {} sits on an occupied site", kind, number));
  occupant_[s] = b;
  site_[b] = s;
  unwrapped_[b] = p;
}

void DllSystem::link(Bead a, Bead b) {
  if (bonds_.full(a) || bonds_.full(b)) {
    throw std::logic_error(std::format("bond {}-{} exceeds {} bonds per bead", a + 1, b + 1, kMaxBonds));
  }
  bonds_.link(a, b);
}

void DllSystem::bond(Bead a, Bead b) {
  if (type(a) != BeadType::Unit || type(b) != BeadType::Unit || a == b || bonded(a, b)) {
    throw std::logic_error(std::format("cannot bond beads {} and {}", a + 1, b + 1));
  }
  tally_.count(bonds_.count(a), bonds_.count(b), [&] { return molecules_.same(a, b); });
  link(a, b);
  molecules_.unite(a, b);
}

StructuralCounts DllSystem::structural_counts() const noexcept {
  StructuralCounts counts{};
  for (std::size_t u = 0; u < states_.size(); ++u) {
    const int code = poly::chem::legacy_type_code(states_[u], bonds_.count(first_unit() + static_cast<Bead>(u)));
    if (code >= 0 && code < poly::chem::kNumLegacyTypeCodes) ++counts[static_cast<std::size_t>(code)];
  }
  return counts;
}

void DllSystem::shift_loop(std::span<const Site> loop, std::span<const std::uint8_t> directions) noexcept {
  const std::size_t k = loop.size();
  const Bead last = occupant_[loop[k - 1]];
  for (std::size_t i = k - 1; i > 0; --i) occupant_[loop[i]] = occupant_[loop[i - 1]];
  occupant_[loop[0]] = last;
  for (std::size_t i = 0; i < k; ++i) {
    const Site to = loop[(i + 1) % k];
    const Bead b = occupant_[to];
    site_[b] = to;
    unwrapped_[b] += FccGeometry::kNeighbors[directions[loop[i]]];
  }
}

std::expected<void, std::string> DllSystem::validate() const {
  const auto fail = [](std::string message) { return std::unexpected(std::move(message)); };
  const std::size_t n = lattice_.num_sites();
  if (occupant_.size() != n || site_.size() != n) return fail("bead and site counts differ");
  for (Site s = 0; s < n; ++s) {
    const Bead b = occupant_[s];
    if (b >= n || site_[b] != s) return fail(std::format("site {} and its bead disagree", s));
    if (lattice_.wrap_any(unwrapped_[b]) != lattice_.position(s)) return fail(std::format("bead {} is not where its site is", b + 1));
  }
  for (Bead b = 0; b < n; ++b) {
    if (type(b) == BeadType::Solvent && bonds_.count(b) != 0) return fail(std::format("solvent bead {} is bonded", b + 1));
    for (const Bead p : bonds(b)) {
      if (p >= n || p == b) return fail(std::format("bead {} has an invalid bond", b + 1));
      if (lattice_.minimum_image(unwrapped_[p] - unwrapped_[b]).norm2() != 2) {
        return fail(std::format("bond {}-{} is broken", b + 1, p + 1));
      }
      if (!bonded(p, b)) return fail(std::format("bond {}-{} is one-sided", b + 1, p + 1));
      if (molecule(p) != molecule(b)) return fail(std::format("bond {}-{} joins two molecules", b + 1, p + 1));
    }
  }
  for (std::size_t u = 0; u < states_.size(); ++u) {
    const UnitState s = states_[u];
    if (s.vinyls > poly::chem::max_vinyls(s.species) || s.radicals > 2) {
      return fail(std::format("unit {} has an impossible chemical state", u + 1));
    }
  }
  return {};
}

}  // namespace dll
