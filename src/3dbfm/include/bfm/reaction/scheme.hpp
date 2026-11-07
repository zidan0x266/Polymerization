// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Reaction schemes and their bookkeeping: the shared FRP/ATRP chemistry of src/common
// (poly/chem/scheme.hpp and poly/chem/progress.hpp), also used by the dll module.
#pragma once

#include "bfm/model/unit.hpp"
#include "poly/chem/progress.hpp"
#include "poly/chem/scheme.hpp"

namespace bfm {

using poly::chem::CompetitionRule;
using poly::chem::EventKind;
using poly::chem::make_scheme;
using poly::chem::parse_competition_rule;
using poly::chem::RadicalChainScheme;
using poly::chem::RateConstants;
using poly::chem::Reaction;
using poly::chem::ReactionProgress;
using poly::chem::ReactionScheme;
using poly::chem::to_string;

}  // namespace bfm
