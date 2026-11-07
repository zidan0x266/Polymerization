// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Trajectory output in the LAMMPS dump format (readable by OVITO, VMD and most analysis tools):
// unwrapped positions in units of the bond length, molecule 0 and type 2 for solvent.
#pragma once

#include <ostream>

#include "dll/model/system.hpp"

namespace dll {

/// Append one frame with the chain beads, and the solvent if `with_solvent`.
void write_dump_frame(std::ostream& out, const DllSystem& system, long step, bool with_solvent);

}  // namespace dll
