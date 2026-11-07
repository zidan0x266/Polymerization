# SPDX-License-Identifier: GPL-2.0-only
# fix poly/react for LAMMPS. This directory is licensed under the GNU General Public License,
# version 2 (see LICENSE), unlike the rest of the toolkit, and uses no code from it.
#
# polyreact_core: the reaction files and the rule engine, tested without LAMMPS.
# polyreactplugin: the fix as a LAMMPS plugin, built only when POLYREACT_LAMMPS_SOURCE names the
# src directory of the LAMMPS it is loaded into (built with -D PKG_PLUGIN=on, the same compiler and
# the same LAMMPS_SIZES).

add_library(polyreact_core STATIC
  lib/reaction_file.cpp
  lib/engine.cpp
  lib/presets.cpp)
target_include_directories(polyreact_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(polyreact_core PRIVATE poly_warnings)
set_target_properties(polyreact_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(polyreact_preset apps/polyreact_preset.cpp)
target_link_libraries(polyreact_preset PRIVATE polyreact_core poly_warnings)

set(POLYREACT_LAMMPS_SOURCE "" CACHE PATH "LAMMPS src directory to build the fix poly/react plugin against")
set(POLYREACT_LAMMPS_SIZES "smallbig" CACHE STRING "LAMMPS_SIZES of that LAMMPS (smallbig or bigbig)")
# LAMMPS formats its messages with its own copy of fmt in C++17 and with std::format from C++20 on,
# so the fix must be compiled with the C++ standard of that LAMMPS (CMAKE_CXX_STANDARD there).
set(POLYREACT_LAMMPS_CXX_STANDARD "17" CACHE STRING "C++ standard of that LAMMPS (17, 20, 23 or 26)")
option(POLYREACT_LAMMPS_MPI "That LAMMPS uses MPI (OFF: the MPI stubs)" ON)

if(POLYREACT_LAMMPS_SOURCE)
  if(NOT EXISTS ${POLYREACT_LAMMPS_SOURCE}/lammpsplugin.h)
    message(FATAL_ERROR "POLYREACT_LAMMPS_SOURCE=${POLYREACT_LAMMPS_SOURCE} has no lammpsplugin.h")
  endif()
  add_library(polyreactplugin MODULE plugin/fix_poly_react.cpp plugin/polyreact_plugin.cpp)
  target_include_directories(polyreactplugin SYSTEM PRIVATE ${POLYREACT_LAMMPS_SOURCE})
  string(TOUPPER ${POLYREACT_LAMMPS_SIZES} sizes)
  target_compile_definitions(polyreactplugin PRIVATE LAMMPS_${sizes} OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)
  target_link_libraries(polyreactplugin PRIVATE polyreact_core poly_warnings)
  if(POLYREACT_LAMMPS_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_link_libraries(polyreactplugin PRIVATE MPI::MPI_CXX)
  else()
    target_include_directories(polyreactplugin SYSTEM PRIVATE ${POLYREACT_LAMMPS_SOURCE}/STUBS)
  endif()
  # the LAMMPS symbols come from the executable that loads the plugin
  if(APPLE)
    target_link_options(polyreactplugin PRIVATE "LINKER:-undefined,dynamic_lookup")
  endif()
  set_target_properties(polyreactplugin PROPERTIES PREFIX "" SUFFIX ".so" CXX_STANDARD ${POLYREACT_LAMMPS_CXX_STANDARD})
endif()
