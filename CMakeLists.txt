cmake_minimum_required(VERSION 3.20)
project(polymerization
  VERSION 2.0.0
  DESCRIPTION "Simulation codes for polymerization and polymer-network formation (Ganesan group, UT Austin)"
  LANGUAGES CXX)

# C++23 is the newest published ISO standard (ISO/IEC 14882:2024) and is fully implemented by
# current compilers. Set -DPOLY_CXX_STANDARD=26 to build against the C++26 working draft.
set(POLY_CXX_STANDARD "23" CACHE STRING "C++ standard to compile with (23 or 26)")
set_property(CACHE POLY_CXX_STANDARD PROPERTY STRINGS 23 26)
set(CMAKE_CXX_STANDARD ${POLY_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(POLY_BUILD_TESTS "Build the unit tests and smoke tests" ON)
option(POLY_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(POLY_NATIVE_ARCH "Optimize for the build machine (-march=native)" OFF)
option(POLY_LTO "Link-time optimization in Release and RelWithDebInfo builds" ON)
set(POLY_PGO "OFF" CACHE STRING "Profile-guided optimization: OFF, GENERATE or USE")
set_property(CACHE POLY_PGO PROPERTY STRINGS OFF GENERATE USE)
set(POLY_PGO_PROFILE "" CACHE FILEPATH "Profile for POLY_PGO=USE (merged .profdata for Clang, profile directory for GCC)")

# Link-time optimization lets the compiler inline across source files (lattice, system, scheme).
if(POLY_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT poly_ipo_supported OUTPUT poly_ipo_message LANGUAGES CXX)
  if(poly_ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
  else()
    message(STATUS "Link-time optimization not available: ${poly_ipo_message}")
  endif()
endif()

# Profile-guided optimization: build with GENERATE, run a representative input, merge the
# profile (Clang: llvm-profdata merge -o run.profdata default.profraw), rebuild with USE.
if(POLY_PGO STREQUAL "GENERATE")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-fprofile-instr-generate)
    add_link_options(-fprofile-instr-generate)
  else()
    add_compile_options(-fprofile-generate)
    add_link_options(-fprofile-generate)
  endif()
elseif(POLY_PGO STREQUAL "USE")
  if(NOT POLY_PGO_PROFILE)
    message(FATAL_ERROR "POLY_PGO=USE needs POLY_PGO_PROFILE")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-fprofile-instr-use=${POLY_PGO_PROFILE})
  else()
    add_compile_options(-fprofile-use=${POLY_PGO_PROFILE} -Wno-missing-profile)
  endif()
elseif(NOT POLY_PGO STREQUAL "OFF")
  message(FATAL_ERROR "POLY_PGO must be OFF, GENERATE or USE")
endif()

# Shared warning/optimization settings for every module.
add_library(poly_warnings INTERFACE)
if(MSVC)
  target_compile_options(poly_warnings INTERFACE /W4 /permissive-)
  if(POLY_WARNINGS_AS_ERRORS)
    target_compile_options(poly_warnings INTERFACE /WX)
  endif()
else()
  # clang folds -Wsign-conversion into -Wconversion; lattice indexing with int coordinates
  # would need a cast on every access, so that one sub-warning is switched off.
  target_compile_options(poly_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
  if(POLY_WARNINGS_AS_ERRORS)
    target_compile_options(poly_warnings INTERFACE -Werror)
  endif()
  if(POLY_NATIVE_ARCH)
    target_compile_options(poly_warnings INTERFACE -march=native)
  endif()
endif()

# Modules: one directory per model under src/.
add_subdirectory(src/common)
add_subdirectory(src/3dbfm)
add_subdirectory(src/randomwalk)
add_subdirectory(src/dll)
# add_subdirectory(src/statistical)   # statistical model, to come
# src/lammps: reactive CG fix, built inside LAMMPS, to come

if(POLY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
