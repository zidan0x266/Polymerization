// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// A set of unit indices with O(1) insertion, removal and random access (shared, src/common).
#pragma once

#include "bfm/core/types.hpp"
#include "poly/index_set.hpp"

namespace bfm {

using IndexSet = poly::IndexSet<Index>;

}  // namespace bfm
