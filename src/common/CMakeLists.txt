# Code shared by every module.
#   poly::common     header-only: random number generator, configuration-file reader, index set,
#                    geometry (lattice vector, FCC lattice), topology (molecules, bond lists) and
#                    output (files, PDB and PSF writers); see README.md
#   poly::chemistry  radical polymerization chemistry (unit states, reaction schemes, bookkeeping
#                    and per-sample tables), used by the reactive lattice modules

add_library(poly_common INTERFACE)
target_include_directories(poly_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(poly::common ALIAS poly_common)

add_library(poly_chemistry STATIC
  lib/scheme.cpp
  lib/tables.cpp)
target_link_libraries(poly_chemistry PUBLIC poly::common PRIVATE poly_warnings)
add_library(poly::chemistry ALIAS poly_chemistry)

install(DIRECTORY include/ DESTINATION include)
