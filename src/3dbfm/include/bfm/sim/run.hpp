// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// The driver: initial placement, athermal relaxation, equilibration, then the diffusion-
// reaction loop of the reactive Monte Carlo scheme. One Monte Carlo step (MCS) consists of N
// (move attempt, reaction attempt) pairs on randomly selected units, exactly as in the legacy
// programs; snapshots are written whenever the vinyl conversion crosses a scheduled level.
// Reaction attempts that cannot involve a reactive unit are skipped without changing the
// probability of any other attempt (see map_reaction_draw).
#pragma once

#include <ostream>

#include "bfm/io/config.hpp"

namespace bfm {

/// Run all samples of `config`, writing the output files and logging to `log`.
void run(Config config, std::ostream& log);

}  // namespace bfm
