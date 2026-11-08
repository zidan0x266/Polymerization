// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/model/initial.hpp"

#include <cstdint>
#include <format>
#include <stdexcept>
#include <utility>

#include "rw/sim/packer.hpp"

namespace dll {

std::string_view to_string(InitialConfiguration c) noexcept {
  switch (c) {
    case InitialConfiguration::RandomWalk: return "random_walk";
    case InitialConfiguration::Aligned: return "aligned";
  }
  return "?";
}

std::optional<InitialConfiguration> parse_initial_configuration(std::string_view s) noexcept {
  if (s == "random_walk") return InitialConfiguration::RandomWalk;
  if (s == "aligned") return InitialConfiguration::Aligned;
  return std::nullopt;
}

std::size_t max_aligned_chains(int cells, std::size_t chain_length) noexcept {
  // The rows along (1, 1, 0): one per even x0 + z, 2n^2 rows of 2n sites each.
  const auto extent = static_cast<std::size_t>(2 * cells);
  return chain_length ? extent * extent / 2 * (extent / chain_length) : 0;
}

namespace {

std::vector<poly::Vec3i> aligned_chains(int cells, std::size_t chain_length, std::size_t n_chains) {
  if (n_chains > max_aligned_chains(cells, chain_length)) {
    throw std::invalid_argument(std::format("at most {} chains of {} beads fit straight into a box of {} cells",
                                            max_aligned_chains(cells, chain_length), chain_length, cells));
  }
  const int extent = 2 * cells;
  const auto per_row = static_cast<std::size_t>(extent) / chain_length;
  std::vector<poly::Vec3i> beads;
  beads.reserve(n_chains * chain_length);
  for (int z = 0; z < extent && beads.size() < n_chains * chain_length; ++z) {
    for (int x0 = z & 1; x0 < extent && beads.size() < n_chains * chain_length; x0 += 2) {
      for (std::size_t c = 0; c < per_row && beads.size() < n_chains * chain_length; ++c) {
        for (std::size_t v = 0; v < chain_length; ++v) {
          const int k = static_cast<int>(c * chain_length + v);
          beads.push_back({x0 + k, k, z});
        }
      }
    }
  }
  return beads;
}

}  // namespace

std::vector<poly::Vec3i> initial_chains(InitialConfiguration kind, int cells, std::size_t chain_length,
                                        std::size_t n_chains, poly::Rng& rng) {
  if (n_chains == 0) return {};
  if (kind == InitialConfiguration::Aligned) return aligned_chains(cells, chain_length, n_chains);
  rw::Packer packer(cells, rw::WalkRules{{60, 90, 120, 180}, rw::SiteWeighting::Uniform, 1000}, rng);
  try {
    return packer.pack(n_chains, chain_length, 0, 1.0).beads;
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(std::format("{}; use initial = aligned for dense systems", e.what()));
  }
}

std::vector<poly::Vec3i> random_free_points(int cells, std::span<const poly::Vec3i> occupied, std::size_t count,
                                             poly::Rng& rng) {
  const poly::FccGeometry lattice(cells);
  std::vector<std::uint8_t> taken(lattice.num_sites(), 0);
  for (const poly::Vec3i p : occupied) taken[lattice.site(lattice.wrap_any(p))] = 1;
  std::vector<poly::Site> free;
  free.reserve(lattice.num_sites() - occupied.size());
  for (poly::Site s = 0; s < lattice.num_sites(); ++s) {
    if (!taken[s]) free.push_back(s);
  }
  if (count > free.size()) {
    throw std::invalid_argument(std::format("{} units do not fit into the {} free sites", count, free.size()));
  }
  // The first `count` entries of a partial Fisher-Yates shuffle are a uniform random sample.
  std::vector<poly::Vec3i> points(count);
  for (std::size_t i = 0; i < count; ++i) {
    const std::size_t j = i + rng.below(static_cast<std::uint32_t>(free.size() - i));
    std::swap(free[i], free[j]);
    points[i] = lattice.position(free[i]);
  }
  return points;
}

}  // namespace dll
