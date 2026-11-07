// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/model/packing.hpp"

#include <format>

namespace rw {

std::expected<void, std::string> validate(const Packing& packing, AngleSet angles) {
  const auto fail = [](std::string message) { return std::unexpected(std::move(message)); };
  if (packing.chain_length < 2) return fail("chains need at least two beads");
  if (packing.size() != packing.chain_beads() + packing.n_crosslinkers) return fail("bead count does not match the composition");

  FccLattice lattice(packing.cells);
  const int extent = lattice.extent();
  const auto modulo = [extent](int c) { return ((c % extent) + extent) % extent; };
  for (std::size_t i = 0; i < packing.size(); ++i) {
    const Vec3i p = packing.beads[i];
    if (((p.x + p.y + p.z) & 1) != 0) return fail(std::format("bead {} is not on a lattice site", i + 1));
    if (!packing.is_chain_bead(i) && (p.x < 0 || p.y < 0 || p.z < 0 || p.x >= extent || p.y >= extent || p.z >= extent)) {
      return fail(std::format("cross-linker {} lies outside the box", i + 1));
    }
    const Site s = lattice.site({modulo(p.x), modulo(p.y), modulo(p.z)});
    if (lattice.occupied(s)) return fail(std::format("bead {} sits on an occupied site", i + 1));
    lattice.occupy(s);
  }

  for (std::size_t c = 0; c < packing.n_chains; ++c) {
    const std::size_t first = c * packing.chain_length;
    for (std::size_t v = 1; v < packing.chain_length; ++v) {
      const Vec3i step = packing.beads[first + v] - packing.beads[first + v - 1];
      if (step.norm2() != 2) return fail(std::format("bond {}-{} is not a nearest-neighbor step", first + v, first + v + 1));
      if (v >= 2 && !angles.allows_steps(packing.beads[first + v - 1] - packing.beads[first + v - 2], step)) {
        return fail(std::format("bond angle at bead {} is not allowed", first + v));
      }
    }
  }
  return {};
}

}  // namespace rw
