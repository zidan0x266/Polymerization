// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Reaction schemes. A scheme is a strategy object: given two neighboring units (or one lone
// unit) it decides stochastically which elementary reaction, if any, takes place. The two
// radical copolymerization mechanisms of the Macromolecules 2025 paper are instances of one
// generic rule set that differs only in which elementary steps exist:
//
//   FRP  : slow initiation, propagation, termination by recombination
//   ATRP : no initiation (all chains live from t=0), propagation, no termination
//
// New mechanisms are added by deriving from ReactionScheme (or by configuring
// RadicalChainScheme::Features).
#pragma once

#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "bfm/core/random.hpp"
#include "bfm/reaction/event.hpp"

namespace bfm {

/// Elementary-step probabilities (per sampled neighbor pair, i.e. per Monte Carlo attempt).
struct RateConstants {
  double k_init{0.0};  ///< initiator decomposition
  double k_p1{0.0};    ///< propagation onto a vinyl of an unreacted monomer / cross-linker
  double k_p2{0.0};    ///< propagation onto a pendant vinyl of a partially reacted cross-linker
  double k_term{0.0};  ///< termination by recombination
};

/// How to resolve a pair for which both propagation and termination are possible (one unit
/// carries a radical, the other both a radical and a pendant vinyl).
enum class CompetitionRule {
  /// Reproduce the legacy implementation bit for bit in distribution: with p = k_p2/(k_p2+k_t)
  /// (or k_p2 if k_p2+k_t < 1) propagate; otherwise terminate with probability
  /// (1-p)*min(1, k_t+p). This is what the published FRP results were generated with.
  Legacy,
  /// A single uniform draw partitions [0,1): [0,k_p2) propagation, [k_p2,k_p2+k_t) termination.
  /// If k_p2 + k_t > 1 the two are renormalized so that one of them always happens.
  Partition
};

[[nodiscard]] std::string_view to_string(CompetitionRule r) noexcept;
[[nodiscard]] std::optional<CompetitionRule> parse_competition_rule(std::string_view s) noexcept;

class ReactionScheme {
 public:
  virtual ~ReactionScheme() = default;

  [[nodiscard]] virtual std::string_view name() const noexcept = 0;

  /// Decide the reaction between a center unit `a` and the neighbor `b` found at a reaction
  /// offset. Returns std::nullopt if nothing happens in this attempt.
  [[nodiscard]] virtual std::optional<ReactionEvent> propose_pair(Index a, UnitView ua, Index b, UnitView ub,
                                                                  Rng& rng) const = 0;

  /// Decide a unimolecular reaction of unit `a` (e.g. initiator decomposition).
  [[nodiscard]] virtual std::optional<ReactionEvent> propose_unimolecular(Index /*a*/, UnitView /*ua*/,
                                                                           Rng& /*rng*/) const {
    return std::nullopt;
  }

  /// Labels of the one-based reaction categories reported by this scheme.
  [[nodiscard]] virtual std::span<const std::string> category_labels() const noexcept = 0;

  /// Initial chemical composition for `n_total` units, in legacy order: initiator-derived
  /// units first, then monovinyl monomers, then divinyl cross-linkers.
  [[nodiscard]] virtual std::vector<UnitState> initial_composition(std::size_t n_total, double monovinyl_fraction,
                                                                   double radical_fraction) const = 0;

  /// Whether the scheme has a unimolecular step at all (lets the driver skip the call).
  [[nodiscard]] virtual bool has_unimolecular_step() const noexcept { return false; }

  /// Whether `u` can drive a reaction: every pair reaction must involve at least one reactive
  /// unit, and only reactive units may undergo a unimolecular step. The driver attempts
  /// reactions only around reactive units (see Simulation), which is exact as long as this
  /// contract holds. The default treats every unit as reactive (always exact, never faster).
  [[nodiscard]] virtual bool is_reactive(UnitView /*u*/) const noexcept { return true; }

  [[nodiscard]] virtual const RateConstants& rate_constants() const noexcept = 0;
};

class RadicalChainScheme final : public ReactionScheme {
 public:
  struct Features {
    bool explicit_initiation{false};  ///< units start as dormant initiators that decompose with k_init
    bool termination{false};          ///< radicals may recombine with k_term
  };

  RadicalChainScheme(std::string name, RateConstants k, Features features, CompetitionRule rule);

  [[nodiscard]] std::string_view name() const noexcept override { return name_; }
  [[nodiscard]] std::optional<ReactionEvent> propose_pair(Index a, UnitView ua, Index b, UnitView ub,
                                                          Rng& rng) const override;
  [[nodiscard]] std::optional<ReactionEvent> propose_unimolecular(Index a, UnitView ua, Rng& rng) const override;
  [[nodiscard]] std::span<const std::string> category_labels() const noexcept override { return labels_; }
  [[nodiscard]] std::vector<UnitState> initial_composition(std::size_t n_total, double monovinyl_fraction,
                                                           double radical_fraction) const override;
  [[nodiscard]] bool has_unimolecular_step() const noexcept override { return features_.explicit_initiation; }
  /// Units with a radical, plus dormant initiators when initiation is explicit.
  [[nodiscard]] bool is_reactive(UnitView u) const noexcept override {
    return u.state.has_radical() || (features_.explicit_initiation && u.state.species == Species::Initiator && u.nbonds == 0);
  }
  [[nodiscard]] const RateConstants& rate_constants() const noexcept override { return k_; }
  [[nodiscard]] const Features& features() const noexcept { return features_; }
  [[nodiscard]] CompetitionRule competition_rule() const noexcept { return rule_; }

  /// Legacy category bookkeeping, exposed for testing.
  [[nodiscard]] static int propagation_category(UnitView attacker, UnitView target) noexcept;
  [[nodiscard]] static int termination_category(UnitView a, UnitView b) noexcept;
  static constexpr int kNumPropagationCategories = 20;
  static constexpr int kNumTerminationCategories = 15;

 private:
  std::string name_;
  RateConstants k_;
  Features features_;
  CompetitionRule rule_;
  std::vector<std::string> labels_;
};

/// Factory for the built-in mechanisms: "frp" and "atrp".
[[nodiscard]] std::unique_ptr<ReactionScheme> make_scheme(std::string_view mechanism, RateConstants k,
                                                          CompetitionRule rule);

}  // namespace bfm
