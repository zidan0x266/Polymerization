// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/sim/reaction.hpp"

#include <array>
#include <numeric>
#include <utility>

namespace dll {

using poly::chem::EventKind;

std::string_view to_string(ReactionPartner p) noexcept {
  return p == ReactionPartner::RandomNeighbor ? "random_neighbor" : "all_neighbors";
}

std::optional<ReactionPartner> parse_reaction_partner(std::string_view s) noexcept {
  if (s == "random_neighbor") return ReactionPartner::RandomNeighbor;
  if (s == "all_neighbors") return ReactionPartner::AllNeighbors;
  return std::nullopt;
}

Reactor::Reactor(const DllSystem& system, const poly::chem::ReactionScheme& scheme, ReactionPartner partner, bool cyclization)
    : scheme_(&scheme), partner_(partner), cyclization_(cyclization), reactive_(system.n_units()), stage_of_(system.n_units(), 0) {
  for (std::size_t u = 0; u < system.n_units(); ++u) update(system, system.first_unit() + static_cast<Bead>(u));
}

void Reactor::update(const DllSystem& system, Bead unit) {
  reactive_.set(unit - system.first_unit(), scheme_->is_reactive(system.view(unit)));
}

bool Reactor::in_stage(const DllSystem& system, Bead b) const noexcept {
  return stage_of_[b - system.first_unit()] == stage_;
}

Bead Reactor::neighbor(const DllSystem& system, Bead b, int k) const noexcept {
  const FccGeometry& lattice = system.lattice();
  return system.occupant(lattice.neighbor(lattice.position(system.site(b)), k));
}

bool Reactor::attempt_pair(DllSystem& system, Bead center, Bead partner, poly::Rng& rng,
                           poly::chem::ReactionProgress& progress) {
  if (system.type(partner) != BeadType::Unit || system.bonded(center, partner)) return false;
  if (!cyclization_ && !system.bonds(center).empty() && !system.bonds(partner).empty() &&
      system.molecule(center) == system.molecule(partner)) {
    return false;
  }
  const auto r = scheme_->propose_pair(system.view(center), system.view(partner), rng);
  if (!r) return false;
  poly::chem::apply(*r, system.state(center), system.state(partner), progress);
  system.bond(center, partner);
  update(system, center);
  update(system, partner);
  ++reactions_;
  return true;
}

bool Reactor::attempt_initiation(DllSystem& system, Bead unit, poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  if (!scheme_->has_unimolecular_step()) return false;
  const auto r = scheme_->propose_unimolecular(system.view(unit), rng);
  if (!r) return false;
  poly::chem::apply(*r, system.state(unit), progress);
  update(system, unit);
  ++reactions_;
  return true;
}

std::size_t Reactor::step(DllSystem& system, poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  reactions_ = 0;
  if (reactive_.empty()) return 0;
  // The units reactive at the start of the stage take part, in random order.
  ++stage_;
  order_.assign(reactive_.members().begin(), reactive_.members().end());
  for (std::size_t i = order_.size(); i > 1; --i) {
    std::swap(order_[i - 1], order_[rng.below(static_cast<std::uint32_t>(i))]);
  }
  for (const std::uint32_t u : order_) stage_of_[u] = stage_;
  if (partner_ == ReactionPartner::RandomNeighbor) {
    step_random_neighbor(system, rng, progress);
  } else {
    step_all_neighbors(system, rng, progress);
  }
  return reactions_;
}

void Reactor::step_random_neighbor(DllSystem& system, poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  // Every unit is the center of one attempt with a random neighbor. An attempt by a unit outside
  // the stage can only succeed toward a unit in the stage, so each unit A in the stage makes two
  // draws: its own attempt, and a random neighbor B that attempts toward A, taken only if B is a
  // unit outside the stage and not reactive (otherwise B's own attempt covers the pair, or B was
  // activated during the stage and waits). Each pair is then attempted with the same
  // probability as if every unit made its attempt.
  for (const std::uint32_t u : order_) {
    const Bead a = system.first_unit() + u;
    const Bead forward = neighbor(system, a, static_cast<int>(rng.below(FccGeometry::kCoordination)));
    if (!attempt_pair(system, a, forward, rng, progress)) attempt_initiation(system, a, rng, progress);
    const Bead backward = neighbor(system, a, static_cast<int>(rng.below(FccGeometry::kCoordination)));
    if (system.type(backward) == BeadType::Unit && !in_stage(system, backward) &&
        !reactive_.contains(backward - system.first_unit())) {
      attempt_pair(system, backward, a, rng, progress);
    }
  }
}

void Reactor::step_all_neighbors(DllSystem& system, poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  // Each unit in the stage tries its neighbors in random order; a pair of two units in the stage
  // is tried only from the one with the lower index, so every pair is tried once.
  std::array<int, FccGeometry::kCoordination> ks{};
  for (const std::uint32_t u : order_) {
    const Bead a = system.first_unit() + u;
    std::iota(ks.begin(), ks.end(), 0);
    for (std::uint32_t i = FccGeometry::kCoordination; i > 1; --i) std::swap(ks[i - 1], ks[rng.below(i)]);
    for (const int k : ks) {
      const Bead b = neighbor(system, a, k);
      if (system.type(b) == BeadType::Unit && in_stage(system, b) && b < a) continue;
      attempt_pair(system, a, b, rng, progress);
    }
    attempt_initiation(system, a, rng, progress);
  }
}

}  // namespace dll
