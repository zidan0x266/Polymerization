// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Starting configurations: precursor chains, and reactive units scattered at random over the
// free sites; the remaining sites are filled with solvent.
#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "poly/fcc.hpp"
#include "poly/random.hpp"

namespace dll {

enum class InitialConfiguration {
  /// Self-avoiding random walks with every bond angle allowed (randomwalk module); works up to
  /// a polymer fraction of about one half.
  RandomWalk,
  /// Straight chains along lattice rows; works up to a fully filled lattice (a melt) but needs
  /// a longer relaxation.
  Aligned
};

[[nodiscard]] std::string_view to_string(InitialConfiguration c) noexcept;
[[nodiscard]] std::optional<InitialConfiguration> parse_initial_configuration(std::string_view s) noexcept;

/// Unwrapped lattice points of `n_chains` chains of `chain_length` beads, chain after chain.
[[nodiscard]] std::vector<poly::Vec3i> initial_chains(InitialConfiguration kind, int cells, std::size_t chain_length,
                                                      std::size_t n_chains, poly::Rng& rng);

/// Most chains of `chain_length` beads that fit straight into lattice rows.
[[nodiscard]] std::size_t max_aligned_chains(int cells, std::size_t chain_length) noexcept;

/// `count` distinct lattice points drawn uniformly from the sites not taken by `occupied`.
[[nodiscard]] std::vector<poly::Vec3i> random_free_points(int cells, std::span<const poly::Vec3i> occupied, std::size_t count,
                                                          poly::Rng& rng);

}  // namespace dll
