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
    : lattice_(cells), chain_length_(chain_length), n_chains_(chain_length ? chains.size() / chain_length : 0) {
  if (chain_length < 2 && !chains.empty()) throw std::invalid_argument("chains need at least two beads");
  if (chain_length && chains.size() % chain_length != 0) throw std::invalid_argument("bead count is not a multiple of the chain length");
  const std::size_t n = lattice_.num_sites();
  if (chains.size() > n) throw std::invalid_argument("more chain beads than lattice sites");

  occupant_.assign(n, kNoBead);
  site_.resize(n);
  unwrapped_.resize(n);
  bonds_.assign(n, {});
  nbonds_.assign(n, 0);

  for (Bead b = 0; b < chains.size(); ++b) {
    const Vec3i p = chains[b];
    if (((p.x + p.y + p.z) & 1) != 0) throw std::invalid_argument(std::format("chain bead {} is not on a lattice site", b + 1));
    const Site s = lattice_.site(lattice_.wrap_any(p));
    if (occupant_[s] != kNoBead) throw std::invalid_argument(std::format("chain bead {} sits on an occupied site", b + 1));
    occupant_[s] = b;
    site_[b] = s;
    unwrapped_[b] = p;
    if (b % chain_length_ != 0) {
      if ((p - chains[b - 1]).norm2() != 2) {
        throw std::invalid_argument(std::format("bond {}-{} is not a nearest-neighbor step", b, b + 1));
      }
      add_bond(b - 1, b);
    }
  }
  auto solvent = static_cast<Bead>(chains.size());
  for (Site s = 0; s < n; ++s) {
    if (occupant_[s] != kNoBead) continue;
    occupant_[s] = solvent;
    site_[solvent] = s;
    unwrapped_[solvent] = lattice_.position(s);
    ++solvent;
  }
}

void DllSystem::add_bond(Bead a, Bead b) noexcept {
  bonds_[a][nbonds_[a]++] = b;
  bonds_[b][nbonds_[b]++] = a;
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
    for (const Bead p : bonds(b)) {
      if (lattice_.minimum_image(unwrapped_[p] - unwrapped_[b]).norm2() != 2) {
        return fail(std::format("bond {}-{} is broken", b + 1, p + 1));
      }
      bool reciprocal = false;
      for (const Bead q : bonds(p)) reciprocal = reciprocal || q == b;
      if (!reciprocal) return fail(std::format("bond {}-{} is one-sided", b + 1, p + 1));
    }
  }
  return {};
}

}  // namespace dll
