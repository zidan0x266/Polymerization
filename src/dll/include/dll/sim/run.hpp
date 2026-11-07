// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Driver: resolve the seed, build or load the starting configuration, equilibrate, then run DLL
// time steps (a motion stage and a reaction stage each) until the conversion reaches
// stop_conversion, nothing can react any more, or the step limit is reached. Snapshots (tables,
// PDB/PSF) are written whenever the conversion passes a scheduled level, as in the 3dbfm
// module. A system without reactive units runs the motion alone for `steps` time steps.
#pragma once

#include <ostream>

#include "dll/io/config.hpp"
#include "dll/model/system.hpp"

namespace dll {

/// Run every sample of `config` and return the system at the end of the last one.
DllSystem run(Config config, std::ostream& log);

}  // namespace dll
