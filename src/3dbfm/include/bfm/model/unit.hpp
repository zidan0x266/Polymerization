// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Chemical state of a lattice unit: the shared definitions of src/common (poly/chem/unit.hpp),
// which also document the legacy "montype" codes.
#pragma once

#include "poly/chem/unit.hpp"

namespace bfm {

using poly::chem::kNumLegacyTypeCodes;
using poly::chem::legacy_residue_name;
using poly::chem::legacy_species_code;
using poly::chem::legacy_type_code;
using poly::chem::max_vinyls;
using poly::chem::Species;
using poly::chem::UnitState;
using poly::chem::UnitView;

}  // namespace bfm
