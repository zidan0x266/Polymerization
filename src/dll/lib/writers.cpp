// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/io/writers.hpp"

#include <cmath>
#include <format>
#include <iterator>
#include <string>

namespace dll {

void write_dump_frame(std::ostream& stream, const DllSystem& system, long step, bool with_solvent) {
  const std::size_t n = with_solvent ? system.size() : system.chain_beads();
  const double unit = 1.0 / std::sqrt(2.0);  // bond length 1
  const double box = unit * system.lattice().extent();
  std::string text;
  text.reserve(64 * n + 256);
  auto out = std::back_inserter(text);
  std::format_to(out, "ITEM: TIMESTEP\n{}\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n", step, n);
  for (int axis = 0; axis < 3; ++axis) std::format_to(out, "0 {:.6f}\n", box);
  std::format_to(out, "ITEM: ATOMS id mol type xu yu zu\n");
  for (Bead b = 0; b < n; ++b) {
    const bool chain = system.type(b) == BeadType::Chain;
    const Vec3i p = system.unwrapped(b);
    std::format_to(out, "{} {} {} {:.6f} {:.6f} {:.6f}\n", b + 1, chain ? b / system.chain_length() + 1 : 0, chain ? 1 : 2,
                   unit * p.x, unit * p.y, unit * p.z);
  }
  stream.write(text.data(), static_cast<std::streamsize>(text.size()));
}

}  // namespace dll
