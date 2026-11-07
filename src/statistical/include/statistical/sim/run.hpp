// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Driver: resolve the seed, build the initial composition, then draw events until the conversion
// reaches stop_conversion or nothing can react any more. Snapshots (tables, PSF) are written
// whenever the conversion passes a scheduled level, as in the lattice modules.
#pragma once

#include <ostream>

#include "statistical/io/config.hpp"
#include "statistical/model/network.hpp"

namespace statistical {

/// Run every sample of `config` and return the network at the end of the last one.
Network run(Config config, std::ostream& log);

}  // namespace statistical
