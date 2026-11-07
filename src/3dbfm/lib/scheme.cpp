// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/reaction/scheme.hpp"

#include <algorithm>
#include <cctype>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace bfm {

std::string_view to_string(CompetitionRule r) noexcept {
  return r == CompetitionRule::Legacy ? "legacy" : "partition";
}

std::optional<CompetitionRule> parse_competition_rule(std::string_view s) noexcept {
  std::string lower(s);
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "legacy") return CompetitionRule::Legacy;
  if (lower == "partition") return CompetitionRule::Partition;
  return std::nullopt;
}

namespace {

const std::vector<std::string> kPropagationLabels = {
    "primary_radical+monomer",          // 1
    "chain_radical+monomer",            // 2
    "crosslinker_radical(4)+monomer",   // 3
    "crosslinker_radical(6)+monomer",   // 4
    "crosslinker_diradical(7)+monomer", // 5
    "primary_radical+crosslinker",      // 6
    "chain_radical+crosslinker",        // 7
    "primary_radical+pendant(4)",       // 8
    "chain_radical+pendant(4)",         // 9
    "primary_radical+pendant(5)",       // 10
    "chain_radical+pendant(5)",         // 11
    "crosslinker_radical(4)+crosslinker",   // 12
    "crosslinker_radical(6)+crosslinker",   // 13
    "crosslinker_diradical(7)+crosslinker", // 14
    "crosslinker_radical(4)+pendant(4)",    // 15
    "crosslinker_radical(4)+pendant(5)",    // 16
    "crosslinker_radical(6)+pendant(4)",    // 17
    "crosslinker_diradical(7)+pendant(4)",  // 18
    "crosslinker_radical(6)+pendant(5)",    // 19
    "crosslinker_diradical(7)+pendant(5)",  // 20
};

const std::vector<std::string> kTerminationLabels = {
    "term:primary+primary",        // 21
    "term:primary+chain",          // 22
    "term:chain+chain",            // 23
    "term:primary+radical(4)",     // 24
    "term:chain+radical(4)",       // 25
    "term:primary+radical(6)",     // 26
    "term:chain+radical(6)",       // 27
    "term:primary+diradical(7)",   // 28
    "term:chain+diradical(7)",     // 29
    "term:radical(4)+radical(4)",  // 30
    "term:radical(4)+radical(6)",  // 31
    "term:radical(4)+diradical(7)",// 32
    "term:radical(6)+radical(6)",  // 33
    "term:radical(6)+diradical(7)",// 34
    "term:diradical(7)+diradical(7)", // 35
};

}  // namespace

RadicalChainScheme::RadicalChainScheme(std::string name, RateConstants k, Features features, CompetitionRule rule)
    : name_(std::move(name)), k_(k), features_(features), rule_(rule) {
  if (k_.k_p1 < 0 || k_.k_p2 < 0 || k_.k_term < 0 || k_.k_init < 0) throw std::invalid_argument("rate constants must be non-negative");
  if (k_.k_p1 > 1 || k_.k_p2 > 1 || k_.k_term > 1 || k_.k_init > 1) throw std::invalid_argument("rate constants are probabilities and must not exceed 1");
  labels_ = kPropagationLabels;
  if (features_.termination) labels_.insert(labels_.end(), kTerminationLabels.begin(), kTerminationLabels.end());
}

std::vector<UnitState> RadicalChainScheme::initial_composition(std::size_t n_total, double monovinyl_fraction,
                                                               double radical_fraction) const {
  if (monovinyl_fraction < 0 || monovinyl_fraction > 1) throw std::invalid_argument("monovinyl_fraction must be in [0,1]");
  if (radical_fraction < 0 || radical_fraction > monovinyl_fraction) {
    throw std::invalid_argument("radical_fraction must be in [0, monovinyl_fraction]");
  }
  // Legacy truncation semantics: N_mono = int(ntot * monovinyl), N_rad = int(ntot * frac_radical).
  const auto n_mono = static_cast<std::size_t>(static_cast<double>(n_total) * monovinyl_fraction);
  const auto n_init = static_cast<std::size_t>(static_cast<double>(n_total) * radical_fraction);
  const UnitState initiator = features_.explicit_initiation ? UnitState::dormant_initiator() : UnitState::primary_radical();
  return std::views::iota(std::size_t{0}, n_total) |
         std::views::transform([&](std::size_t i) {
           if (i < n_init) return initiator;
           return i < n_mono ? UnitState::monovinyl_monomer() : UnitState::divinyl_monomer();
         }) |
         std::ranges::to<std::vector>();
}

int RadicalChainScheme::propagation_category(UnitView attacker, UnitView target) noexcept {
  const int att = attacker.legacy_code();
  const int tgt = target.legacy_code();
  // columns: target code 0 (monomer), 3 (cross-linker), 4 (pendant on radical unit), 5 (pendant)
  const auto column = [tgt]() {
    switch (tgt) {
      case 0: return 0;
      case 3: return 1;
      case 4: return 2;
      case 5: return 3;
      default: return -1;
    }
  }();
  if (column < 0) return 0;
  switch (att) {
    case 1: {
      constexpr int fresh[4] = {1, 6, 8, 10};
      return attacker.nbonds == 0 ? fresh[column] : fresh[column] + 1;
    }
    case 4: {
      constexpr int table[4] = {3, 12, 15, 16};
      return table[column];
    }
    case 6: {
      constexpr int table[4] = {4, 13, 17, 19};
      return table[column];
    }
    case 7: {
      constexpr int table[4] = {5, 14, 18, 20};
      return table[column];
    }
    default: return 0;
  }
}

int RadicalChainScheme::termination_category(UnitView a, UnitView b) noexcept {
  int ca = a.legacy_code(), cb = b.legacy_code();
  std::uint8_t na = a.nbonds, nb = b.nbonds;
  if (ca > cb) {
    std::swap(ca, cb);
    std::swap(na, nb);
  }
  const int base = kNumPropagationCategories;
  if (ca == 1 && cb == 1) {
    if (na == 0 && nb == 0) return base + 1;
    if (na == 0 || nb == 0) return base + 2;
    return base + 3;
  }
  if (ca == 1) {
    // the type-1 unit is `a` after sorting; split on whether it is a fresh primary radical
    const int offset = (na == 0) ? 0 : 1;
    switch (cb) {
      case 4: return base + 4 + offset;
      case 6: return base + 6 + offset;
      case 7: return base + 8 + offset;
      default: return 0;
    }
  }
  if (ca == 4 && cb == 4) return base + 10;
  if (ca == 4 && cb == 6) return base + 11;
  if (ca == 4 && cb == 7) return base + 12;
  if (ca == 6 && cb == 6) return base + 13;
  if (ca == 6 && cb == 7) return base + 14;
  if (ca == 7 && cb == 7) return base + 15;
  return 0;
}

std::optional<ReactionEvent> RadicalChainScheme::propose_pair(Index a, UnitView ua, Index b, UnitView ub,
                                                              Rng& rng) const {
  // Direction of a possible propagation: the center unit attacks first (legacy order).
  Index attacker = kNone, target = kNone;
  UnitView vatt{}, vtgt{};
  if (ua.state.has_radical() && ub.state.has_vinyl()) {
    attacker = a; target = b; vatt = ua; vtgt = ub;
  } else if (ub.state.has_radical() && ua.state.has_vinyl()) {
    attacker = b; target = a; vatt = ub; vtgt = ua;
  }
  const bool can_propagate = attacker != kNone;
  const bool can_terminate = features_.termination && ua.state.has_radical() && ub.state.has_radical();
  if (!can_propagate && !can_terminate) return std::nullopt;

  const auto propagation = [&] {
    return ReactionEvent{EventKind::Propagation, attacker, target, propagation_category(vatt, vtgt)};
  };
  const auto termination = [&] { return ReactionEvent{EventKind::Termination, a, b, termination_category(ua, ub)}; };

  if (can_propagate && !can_terminate) {
    const double p = vtgt.state.is_unreacted_monomer() ? k_.k_p1 : k_.k_p2;
    if (rng.bernoulli(p)) return propagation();
    return std::nullopt;
  }
  if (!can_propagate) {
    if (rng.bernoulli(k_.k_term)) return termination();
    return std::nullopt;
  }

  // Both are possible: the target vinyl is necessarily a pendant (a unit with a radical cannot
  // be an unreacted monomer), so the propagation constant is k_p2.
  switch (rule_) {
    case CompetitionRule::Legacy: {
      const double sum = k_.k_p2 + k_.k_term;
      const double p_prop = (sum >= 1.0) ? k_.k_p2 / sum : k_.k_p2;
      if (rng.uniform01() < p_prop) return propagation();
      const double r2 = rng.uniform01();
      const double r3 = rng.uniform01();
      if (r2 >= p_prop && r3 < k_.k_term + p_prop) return termination();
      return std::nullopt;
    }
    case CompetitionRule::Partition: {
      const double sum = k_.k_p2 + k_.k_term;
      const double scale = (sum > 1.0) ? 1.0 / sum : 1.0;
      const double u = rng.uniform01();
      if (u < k_.k_p2 * scale) return propagation();
      if (u < sum * scale) return termination();
      return std::nullopt;
    }
  }
  std::unreachable();
}

std::optional<ReactionEvent> RadicalChainScheme::propose_unimolecular(Index a, UnitView ua, Rng& rng) const {
  if (!features_.explicit_initiation) return std::nullopt;
  if (ua.state.species != Species::Initiator || ua.state.has_radical() || ua.nbonds != 0) return std::nullopt;
  if (rng.bernoulli(k_.k_init)) return ReactionEvent{EventKind::Initiation, a, kNone, 0};
  return std::nullopt;
}

std::unique_ptr<ReactionScheme> make_scheme(std::string_view mechanism, RateConstants k, CompetitionRule rule) {
  std::string lower(mechanism);
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "frp") {
    return std::make_unique<RadicalChainScheme>("FRP", k, RadicalChainScheme::Features{true, true}, rule);
  }
  if (lower == "atrp") {
    return std::make_unique<RadicalChainScheme>("ATRP", k, RadicalChainScheme::Features{false, false}, rule);
  }
  throw std::invalid_argument("unknown reaction mechanism '" + std::string(mechanism) + "' (expected frp or atrp)");
}

}  // namespace bfm
