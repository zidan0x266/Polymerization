// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "statistical/sim/kinetics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <utility>

namespace statistical {

using poly::chem::EventKind;
using poly::chem::RadicalChainScheme;
using poly::chem::Reaction;
using poly::chem::Species;

namespace {

void erase_at(std::vector<Unit>& v, std::size_t i) noexcept {
  v[i] = v.back();
  v.pop_back();
}

std::uint32_t draw(poly::Rng& rng, std::size_t n) { return rng.below(static_cast<std::uint32_t>(n)); }

}  // namespace

Kinetics::Kinetics(Network& network, RadicalChainScheme::Features features, poly::chem::RateConstants rates, bool cyclization)
    : network_(&network), features_(features), rates_(rates), cyclization_(cyclization), volume_(static_cast<double>(network.size())) {
  for (Unit u = 0; u < network.size(); ++u) {
    const UnitState s = network.state(u);
    for (int r = 0; r < s.radicals; ++r) radicals_.push_back(u);
    if (s.is_unreacted_monomer()) (s.species == Species::Divinyl ? crosslinkers_ : monomers_).push_back(u);
    for (int v = 0; v < pendant_vinyls(s); ++v) pendants_.push_back(u);
    if (features_.explicit_initiation && s.species == Species::Initiator && !s.has_radical() && network.bonds(u).empty()) {
      dormant_.push_back(u);
    }
  }
}

bool Kinetics::allowed(Unit a, Unit b) const noexcept {
  return a != b && !network_->bonded(a, b) && (cyclization_ || network_->molecule(a) != network_->molecule(b));
}

void Kinetics::react(Unit a, Unit b, const Reaction& r, poly::chem::ReactionProgress& progress) {
  UnitState sa = network_->state(a), sb = network_->state(b);
  poly::chem::apply(r, sa, sb, progress);
  network_->set_state(a, sa);
  network_->set_state(b, sb);
  network_->bond(a, b);
}

Event Kinetics::step(poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  const auto r = static_cast<double>(radicals_.size());
  const std::array<double, 5> rate{
      features_.explicit_initiation ? rates_.k_init * static_cast<double>(dormant_.size()) : 0.0,
      rates_.k_p1 * r * static_cast<double>(monomers_.size()) / volume_,
      rates_.k_px.value_or(rates_.k_p1) * r * static_cast<double>(crosslinkers_.size()) / volume_,
      rates_.k_p2 * r * static_cast<double>(pendants_.size()) / volume_,
      features_.termination ? rates_.k_term * r * (r - 1.0) / 2.0 / volume_ : 0.0};
  double total = 0.0;
  for (const double a : rate) total += a;
  if (!(total > 0.0)) return Event::Exhausted;
  time_ += -std::log(1.0 - rng.uniform01()) / total;

  // The event class: the first whose cumulative rate exceeds the draw (rounding: the last one
  // with a positive rate).
  const double x = rng.uniform01() * total;
  int kind = 0;
  double cumulative = 0.0;
  for (int k = 0; k < static_cast<int>(rate.size()); ++k) {
    if (rate[k] <= 0.0) continue;
    kind = k;
    cumulative += rate[k];
    if (x < cumulative) break;
  }
  switch (kind) {
    case 0: return initiate(rng, progress);
    case 1: return propagate(monomers_, rng, progress);
    case 2: return propagate(crosslinkers_, rng, progress);
    case 3: return propagate_pendant(rng, progress);
    default: return terminate(rng, progress);
  }
}

Event Kinetics::initiate(poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  const std::size_t j = draw(rng, dormant_.size());
  const Unit u = dormant_[j];
  erase_at(dormant_, j);
  UnitState s = network_->state(u);
  poly::chem::apply(Reaction{EventKind::Initiation, false, 0}, s, progress);
  network_->set_state(u, s);
  radicals_.push_back(u);
  return Event::Reaction;
}

Event Kinetics::propagate(std::vector<Unit>& targets, poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  // An unreacted monomer or cross-linker is a molecule of its own, so every pair is allowed.
  const std::size_t i = draw(rng, radicals_.size());
  const std::size_t j = draw(rng, targets.size());
  const Unit a = radicals_[i], t = targets[j];
  const Reaction r{EventKind::Propagation, false, RadicalChainScheme::propagation_category(network_->view(a), network_->view(t))};
  erase_at(radicals_, i);
  erase_at(targets, j);
  react(a, t, r, progress);
  radicals_.push_back(t);
  if (pendant_vinyls(network_->state(t)) > 0) pendants_.push_back(t);  // a cross-linker keeps a pendant vinyl
  return Event::Reaction;
}

Event Kinetics::propagate_pendant(poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  const std::size_t i = draw(rng, radicals_.size());
  const Unit a = radicals_[i];
  // Pendant vinyls the radical may not take: its own, those of its bonded neighbors (no double
  // bonds), and, without cyclization, all of its molecule. Partners are drawn until one is allowed.
  std::size_t excluded = 0;
  if (cyclization_) {
    excluded = static_cast<std::size_t>(pendant_vinyls(network_->state(a)));
    for (const Unit b : network_->bonds(a)) excluded += static_cast<std::size_t>(pendant_vinyls(network_->state(b)));
  } else {
    excluded = network_->molecule_pendants(a);
  }
  if (excluded >= pendants_.size()) return Event::Idle;
  std::size_t j = draw(rng, pendants_.size());
  while (!allowed(a, pendants_[j])) j = draw(rng, pendants_.size());
  const Unit t = pendants_[j];
  const Reaction r{EventKind::Propagation, false, RadicalChainScheme::propagation_category(network_->view(a), network_->view(t))};
  erase_at(radicals_, i);
  erase_at(pendants_, j);
  react(a, t, r, progress);
  radicals_.push_back(t);
  return Event::Reaction;
}

Event Kinetics::terminate(poly::Rng& rng, poly::chem::ReactionProgress& progress) {
  const std::size_t i = draw(rng, radicals_.size());
  const Unit a = radicals_[i];
  std::size_t excluded = 0;
  if (cyclization_) {
    excluded = network_->state(a).radicals;
    for (const Unit b : network_->bonds(a)) excluded += network_->state(b).radicals;
  } else {
    excluded = network_->molecule_radicals(a);
  }
  if (excluded >= radicals_.size()) return Event::Idle;
  std::size_t j = draw(rng, radicals_.size());
  while (!allowed(a, radicals_[j])) j = draw(rng, radicals_.size());
  const Unit b = radicals_[j];
  const Reaction r{EventKind::Termination, false, RadicalChainScheme::termination_category(network_->view(a), network_->view(b))};
  erase_at(radicals_, std::max(i, j));
  erase_at(radicals_, std::min(i, j));
  react(a, b, r, progress);
  return Event::Reaction;
}

std::expected<void, std::string> Kinetics::validate() const {
  std::vector<Unit> radicals, monomers, crosslinkers, pendants, dormant;
  for (Unit u = 0; u < network_->size(); ++u) {
    const UnitState s = network_->state(u);
    for (int r = 0; r < s.radicals; ++r) radicals.push_back(u);
    if (s.is_unreacted_monomer()) (s.species == Species::Divinyl ? crosslinkers : monomers).push_back(u);
    for (int v = 0; v < pendant_vinyls(s); ++v) pendants.push_back(u);
    if (features_.explicit_initiation && s.species == Species::Initiator && !s.has_radical() && network_->bonds(u).empty()) {
      dormant.push_back(u);
    }
  }
  const auto same = [](std::vector<Unit> expected, std::vector<Unit> actual) {
    std::ranges::sort(expected);
    std::ranges::sort(actual);
    return expected == actual;
  };
  if (!same(radicals, radicals_)) return std::unexpected(std::string("the radical list is out of date"));
  if (!same(monomers, monomers_) || !same(crosslinkers, crosslinkers_)) return std::unexpected(std::string("the monomer lists are out of date"));
  if (!same(pendants, pendants_)) return std::unexpected(std::string("the pendant list is out of date"));
  if (!same(dormant, dormant_)) return std::unexpected(std::string("the initiator list is out of date"));
  return {};
}

}  // namespace statistical
