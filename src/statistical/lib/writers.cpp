// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "statistical/io/writers.hpp"

#include "poly/chem/psf.hpp"

namespace statistical {

void write_psf(const std::filesystem::path& path, const Network& network) {
  poly::chem::write_psf(
      path, network.size(), [&](std::size_t u) { return poly::chem::legacy_residue_name(network.view(static_cast<Unit>(u)).legacy_code()); },
      [&](std::size_t u) { return network.bonds(static_cast<Unit>(u)); });
}

}  // namespace statistical
