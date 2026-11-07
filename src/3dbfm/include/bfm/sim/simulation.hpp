// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The driver: initial placement, athermal relaxation, equilibration, then the diffusion-
// reaction loop of the reactive Monte Carlo scheme. One Monte Carlo step (MCS) consists of N
// (move attempt, reaction attempt) pairs on randomly selected units, exactly as in the legacy
// programs; snapshots are written whenever the vinyl conversion crosses a scheduled level.
// Reaction attempts that cannot involve a reactive unit are skipped without changing the
// probability of any other attempt (see map_reaction_draw).
#pragma once

#include <chrono>
#include <cstdint>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "bfm/core/index_set.hpp"
#include "bfm/core/random.hpp"
#include "bfm/io/config.hpp"
#include "bfm/lattice/bond_set.hpp"
#include "bfm/model/system.hpp"
#include "bfm/reaction/scheme.hpp"

namespace bfm {

/// Running totals of the reaction phase for one sample.
struct ReactionProgress {
  std::uint64_t total_vinyls{0};     ///< N_mono + 2 N_di + N_initiator (legacy vinyl budget)
  std::uint64_t consumed_vinyls{0};  ///< counts the initiator-derived units as consumed (legacy definition)
  std::uint64_t propagations{0};
  std::uint64_t terminations{0};
  std::uint64_t initiations{0};
  std::uint64_t live_radicals{0};
  std::vector<std::uint64_t> category_counts;  ///< one-based categories; index 0 unused

  [[nodiscard]] double conversion() const noexcept {
    return total_vinyls ? static_cast<double>(consumed_vinyls) / static_cast<double>(total_vinyls) : 0.0;
  }
};

/// A legacy reaction attempt: the center unit and the reaction-offset slot taken from it.
struct ReactionAttempt {
  Index center{kNone};
  int slot{0};
  constexpr bool operator==(const ReactionAttempt&) const noexcept = default;
};

/// Map a uniform draw (x, slot), x in [0, n) and slot in [0, 54), to the legacy attempt
/// (center, slot) it stands for, or to nothing when no reactive unit can be involved. With A the
/// set of reactive units:
///   x <  |A|          center A[x] at `slot` (the legacy attempt from a reactive center);
///   |A| <= x < 2|A|   the reactive unit A[x - |A|] is the partner: the center is the unit at
///                     offset `slot` from it, taken only if it is not reactive itself (otherwise
///                     the forward case covers the pair);
///   x >= 2|A|         nothing.
/// Every legacy attempt that involves a reactive unit is produced by exactly one draw, so its
/// probability stays 1/(54 n). When 2|A| > n the draw is used as in the legacy code.
[[nodiscard]] std::optional<ReactionAttempt> map_reaction_draw(const BfmSystem& system, const IndexSet& reactive,
                                                               std::uint32_t x, int slot) noexcept;

class Simulation {
 public:
  Simulation(Config config, std::unique_ptr<ReactionScheme> scheme, std::ostream& log);

  /// Run all samples configured in `config`.
  void run();

  /// Run one sample (numbered as in the output files).
  void run_sample(int sample);

  [[nodiscard]] const Config& config() const noexcept { return config_; }
  [[nodiscard]] const ReactionScheme& scheme() const noexcept { return *scheme_; }
  [[nodiscard]] const BfmSystem& system() const noexcept { return *system_; }
  [[nodiscard]] const ReactionProgress& progress() const noexcept { return progress_; }

 private:
  void prepare_system(int sample);
  void relax(long n_mcs, const char* label);
  void react(int sample);
  bool attempt_reaction(Index center, int slot);
  void apply(const ReactionEvent& ev, Index center, Index partner, int slot);
  void rebuild_reactive_set();
  void update_reactive(Index id);
  void snapshot(int sample, long mcs, double level);
  void log_line(const std::string& message) const;

  Config config_;
  std::unique_ptr<ReactionScheme> scheme_;
  std::ostream& log_;
  BondVectorSet bonds_;
  std::unique_ptr<BfmSystem> system_;
  std::unique_ptr<Rng> rng_;
  ReactionProgress progress_{};
  IndexSet reactive_;  ///< units for which scheme_->is_reactive() holds
  std::chrono::steady_clock::time_point started_;
};

}  // namespace bfm
