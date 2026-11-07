# SPDX-License-Identifier: GPL-2.0-only
# fix poly/react for LAMMPS. This directory is licensed under the GNU General Public License,
# version 2 (see LICENSE), unlike the rest of the toolkit, and uses no code from it.
# Phase 1: the reaction files and the rule engine (library polyreact_core) with their tests; the
# fix itself follows.

add_library(polyreact_core STATIC
  lib/reaction_file.cpp
  lib/engine.cpp
  lib/presets.cpp)
target_include_directories(polyreact_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(polyreact_core PRIVATE poly_warnings)

add_executable(polyreact_preset apps/polyreact_preset.cpp)
target_link_libraries(polyreact_preset PRIVATE polyreact_core poly_warnings)
