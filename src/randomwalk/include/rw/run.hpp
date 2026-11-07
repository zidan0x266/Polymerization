// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Driver: resolve the seed, pack the chains and cross-linkers, check the result and write the
// LAMMPS data file, the XYZ file and run_config.ini to the output directory.
#pragma once

#include <ostream>

#include "rw/io/config.hpp"
#include "rw/model/packing.hpp"

namespace rw {

Packing run(Config config, std::ostream& log);

}  // namespace rw
