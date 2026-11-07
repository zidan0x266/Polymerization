// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Reaction files for radical copolymerization of monovinyl monomers and divinyl cross-linkers,
// with the chemistry of the 3DBFM module (Z. Zhang et al., Macromolecules 2025, 58, 3168): free
// radical (FRP: initiation, propagation, termination by recombination) and atom transfer radical
// polymerization (ATRP: all chains start at once, no termination). Each bead is a monomer unit;
// its type is its structural unit:
//
//   MON  unreacted monovinyl monomer          FOU  cross-linker with a radical and a pendant vinyl
//   PRI  initiator fragment with its radical  FIV  cross-linker with a pendant vinyl
//   ONE  monovinyl unit with a radical        SIX  cross-linker with a radical, both vinyls reacted
//   TWO  saturated unit or used initiator     SEV  cross-linker with two radicals
//   THR  unreacted divinyl cross-linker       EIG  fully reacted cross-linker
//   INI  dormant initiator (FRP)
//
// The reactions follow the reaction categories of the original code (propagation 1-20, then
// termination 21-35, then initiation), so the counts per reaction can be compared column by column.
#pragma once

#include <array>
#include <optional>
#include <string>
#include <string_view>

#include "polyreact/reaction_file.hpp"

namespace polyreact {

enum class Mechanism { Frp, Atrp };

[[nodiscard]] std::optional<Mechanism> parse_mechanism(std::string_view s) noexcept;

inline constexpr std::array<std::string_view, 11> kRadicalTypes{"MON", "PRI", "ONE", "TWO", "THR", "FOU",
                                                                "FIV", "SIX", "SEV", "EIG", "INI"};

/// Probabilities per check, as the rate constants of the 3DBFM module.
struct RadicalRates {
  double k_init{0.0};  ///< initiator decomposition (FRP)
  double k_p1{0.0};    ///< propagation onto an unreacted monomer
  double k_p2{0.0};    ///< propagation onto a pendant vinyl
  double k_term{0.0};  ///< termination by recombination (FRP)
  std::optional<double> k_px{};  ///< propagation onto an unreacted cross-linker; unset: k_p1
};

struct PresetOptions {
  long interval{1};
  double cutoff{1.2};
  int bond_type{1};
  std::optional<int> angle_type;
  Competition competition{Competition::Legacy};
  bool cyclization{true};  ///< allow bonds within a molecule (rings, intramolecular cross-links)
};

/// The reaction file for `mechanism`; throws std::invalid_argument for rates outside [0, 1].
[[nodiscard]] std::string radical_reaction_file(Mechanism mechanism, const RadicalRates& rates, const PresetOptions& options);

}  // namespace polyreact
