// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Driver: resolve the seed, build the starting configuration, run the DLL time steps and write
// the log, the trajectory and run_config.ini.
#pragma once

#include <ostream>

#include "dll/io/config.hpp"
#include "dll/model/system.hpp"

namespace dll {

/// Run a simulation and return the final system.
DllSystem run(Config config, std::ostream& log);

}  // namespace dll
