// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/model/system.hpp"

#include <numeric>
#include <string>

namespace bfm {

BfmSystem::BfmSystem(PeriodicLattice lattice, const BondVectorSet& bonds)
    : lattice_(std::move(lattice)), bonds_(&bonds) {}

void BfmSystem::place_randomly(std::span<const UnitState> composition, Rng& rng) {
  const auto n = static_cast<Index>(composition.size());
  const std::size_t capacity = static_cast<std::size_t>(lattice_.lx() / 2) * (lattice_.ly() / 2) * (lattice_.lz() / 2);
  if (static_cast<std::size_t>(n) > capacity) {
    throw std::invalid_argument("too many units for the lattice: at most L^3/8 cubes fit on the even sub-lattice");
  }
  lattice_.clear();
  monomers_.assign(static_cast<std::size_t>(n), Monomer{});
  images_.assign(static_cast<std::size_t>(n), Vec3i{});
  parent_.resize(static_cast<std::size_t>(n));
  std::iota(parent_.begin(), parent_.end(), Index{0});
  tally_ = {};

  const auto hx = static_cast<std::uint32_t>(lattice_.lx() / 2);
  const auto hy = static_cast<std::uint32_t>(lattice_.ly() / 2);
  const auto hz = static_cast<std::uint32_t>(lattice_.lz() / 2);
  for (Index id = 0; id < n; ++id) {
    while (true) {
      const Vec3i anchor{2 * static_cast<int>(rng.below(hx)), 2 * static_cast<int>(rng.below(hy)),
                         2 * static_cast<int>(rng.below(hz))};
      if (lattice_.occupied(anchor)) continue;
      lattice_.occupy_cube(anchor, id);
      auto& m = monomers_[static_cast<std::size_t>(id)];
      m.set_anchor(anchor);
      m.state = composition[static_cast<std::size_t>(id)];
      break;
    }
  }
}

void BfmSystem::load(std::vector<Monomer> monomers, std::span<const Vec3i> unwrapped) {
  lattice_.clear();
  monomers_ = std::move(monomers);
  const auto n = size();
  if (!unwrapped.empty() && unwrapped.size() != monomers_.size()) {
    throw std::invalid_argument("restart configuration: unwrapped coordinates do not match the number of units");
  }
  images_.assign(static_cast<std::size_t>(n), Vec3i{});
  parent_.resize(static_cast<std::size_t>(n));
  std::iota(parent_.begin(), parent_.end(), Index{0});
  tally_ = {};
  for (Index id = 0; id < n; ++id) {
    const auto& m = monomers_[static_cast<std::size_t>(id)];
    const Vec3i a = m.anchor();
    if (a.x < 0 || a.y < 0 || a.z < 0 || a.x >= lattice_.lx() || a.y >= lattice_.ly() || a.z >= lattice_.lz()) {
      throw std::runtime_error("restart configuration has a unit outside the box (unit " + std::to_string(id) + ")");
    }
    if (!lattice_.cube_is_free(a)) {
      throw std::runtime_error("restart configuration has overlapping units (unit " + std::to_string(id) + ")");
    }
    lattice_.occupy_cube(a, id);
    if (!unwrapped.empty()) {
      // The unwrapped anchor differs from the wrapped one by whole box lengths.
      const Vec3i d = unwrapped[static_cast<std::size_t>(id)] - a;
      if (d.x % lattice_.lx() != 0 || d.y % lattice_.ly() != 0 || d.z % lattice_.lz() != 0) {
        throw std::runtime_error("restart configuration: unwrapped and wrapped anchors disagree (unit " + std::to_string(id) + ")");
      }
      images_[static_cast<std::size_t>(id)] = {d.x / lattice_.lx(), d.y / lattice_.ly(), d.z / lattice_.lz()};
    }
  }
  for (Index id = 0; id < n; ++id) {
    for (const Bond b : monomers_[static_cast<std::size_t>(id)].bonds()) unite(id, b.partner);
  }
  validate();
}

bool BfmSystem::attempt_move(Index id, Direction d) noexcept {
  Monomer& m = monomers_[static_cast<std::size_t>(id)];
  std::array<BondIndex, kMaxBondsPerUnit> shifted{};
  for (int i = 0; i < m.nbonds; ++i) {
    const auto s = bonds_->shifted(m.vectors[i], d);
    if (!s) return false;  // the bond would leave the allowed set
    shifted[i] = *s;
  }
  const Vec3i from = m.anchor();
  if (!lattice_.front_is_free(from, d)) return false;  // excluded volume

  for (int i = 0; i < m.nbonds; ++i) {
    Monomer& p = monomers_[static_cast<std::size_t>(m.partners[i])];
    p.vectors[p.slot_of(id)] = BondVectorSet::reverse(shifted[i]);
    m.vectors[i] = shifted[i];
  }
  const Vec3i to = lattice_.shift_cube(from, d, id);
  m.set_anchor(to);
  const Vec3i step = unit_step(d);
  if (to - from != step) images_[static_cast<std::size_t>(id)] += step;  // crossed the periodic boundary
  return true;
}

void BfmSystem::add_bond(Index a, Index b, BondIndex vector_from_b_to_a) {
  Monomer& ma = monomers_[static_cast<std::size_t>(a)];
  Monomer& mb = monomers_[static_cast<std::size_t>(b)];
  if (ma.nbonds >= kMaxBondsPerUnit || mb.nbonds >= kMaxBondsPerUnit) {
    throw std::logic_error("attempt to exceed the maximum number of bonds per unit");
  }
  ma.add_bond(Bond{b, vector_from_b_to_a});
  mb.add_bond(Bond{a, BondVectorSet::reverse(vector_from_b_to_a)});
}

void BfmSystem::bond_through_offset(Index center, Index partner, int slot) {
  const Monomer& mc = monomers_[static_cast<std::size_t>(center)];
  const Monomer& mp = monomers_[static_cast<std::size_t>(partner)];
  if (mc.nbonds == 0 && mp.nbonds == 0) {
    ++tally_.monomer_monomer;
  } else if (mc.nbonds != 0 && mp.nbonds != 0) {
    if (same_molecule(center, partner)) {
      ++tally_.intramolecular;
    } else {
      ++tally_.polymer_polymer;
    }
  } else {
    ++tally_.polymer_monomer;
  }
  // The partner sits at center + offset, so the vector from partner to center is -offset.
  add_bond(center, partner, bonds_->reaction_offsets()[slot].bond_from_partner_to_center);
  unite(center, partner);
}

Index BfmSystem::molecule_root(Index id) noexcept {
  while (parent_[static_cast<std::size_t>(id)] != id) {
    const Index grand = parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(id)])];
    parent_[static_cast<std::size_t>(id)] = grand;  // path halving
    id = grand;
  }
  return id;
}

void BfmSystem::unite(Index a, Index b) noexcept {
  const Index ra = molecule_root(a);
  const Index rb = molecule_root(b);
  if (ra != rb) parent_[static_cast<std::size_t>(rb)] = ra;
}

std::array<std::int64_t, kNumLegacyTypeCodes> BfmSystem::structural_counts() const noexcept {
  std::array<std::int64_t, kNumLegacyTypeCodes> counts{};
  for (const Monomer& m : monomers_) {
    const int code = legacy_type_code(m.state, m.nbonds);
    if (code >= 0 && code < kNumLegacyTypeCodes) ++counts[static_cast<std::size_t>(code)];
  }
  return counts;
}

void BfmSystem::validate() const {
  const Index n = size();
  // Every site of every cube must be occupied, no site may be claimed twice, and the lattice
  // must hold no occupied site or anchor that belongs to no unit.
  std::vector<std::uint8_t> claimed(lattice_.site_count(), 0);
  for (Index id = 0; id < n; ++id) {
    const Monomer& m = monomers_[static_cast<std::size_t>(id)];
    const Vec3i a = m.anchor();
    if (lattice_.anchored_unit(a) != id) throw std::runtime_error("lattice anchor mismatch at unit " + std::to_string(id));
    for (int dx = 0; dx < 2; ++dx) {
      for (int dy = 0; dy < 2; ++dy) {
        for (int dz = 0; dz < 2; ++dz) {
          const Vec3i p = lattice_.wrap(a + Vec3i{dx, dy, dz});
          if (!lattice_.occupied(p) || claimed[lattice_.site_index(p)]++ != 0) {
            throw std::runtime_error("lattice occupancy mismatch at unit " + std::to_string(id));
          }
        }
      }
    }
  }
  if (lattice_.occupied_count() != 8 * static_cast<std::size_t>(n) || lattice_.anchor_count() != static_cast<std::size_t>(n)) {
    throw std::runtime_error("lattice holds sites or anchors that belong to no unit");
  }

  for (Index id = 0; id < n; ++id) {
    const Monomer& m = monomers_[static_cast<std::size_t>(id)];
    for (const Bond b : m.bonds()) {
      if (b.partner < 0 || b.partner >= n || b.partner == id) throw std::runtime_error("invalid bond partner");
      const Monomer& p = monomers_[static_cast<std::size_t>(b.partner)];
      const int back = p.slot_of(id);
      if (back < 0 || p.vectors[back] != BondVectorSet::reverse(b.vector)) {
        throw std::runtime_error("non-reciprocal bond between units " + std::to_string(id) + " and " + std::to_string(b.partner));
      }
      // vector(b) must equal the minimum-image displacement from the partner to this unit
      Vec3i disp = m.anchor() - p.anchor();
      for (int ax = 0; ax < 3; ++ax) {
        const int L = lattice_.extent(ax);
        int& c = (ax == 0 ? disp.x : (ax == 1 ? disp.y : disp.z));
        if (c > L / 2) c -= L;
        if (c < -L / 2) c += L;
      }
      if (disp != bonds_->vector(b.vector)) {
        throw std::runtime_error("bond vector does not match geometry for units " + std::to_string(id) + " and " + std::to_string(b.partner));
      }
    }
  }
}

}  // namespace bfm
