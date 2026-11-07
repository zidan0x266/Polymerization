# randomwalk module: packs precursor chains (self-avoiding random walks on an FCC lattice) and
# cross-linkers into a box as the starting point of end-linking simulations in LAMMPS. Builds
# the static library `randomwalk` and the command-line front end `rw_pack`.

add_library(randomwalk STATIC
  lib/fcc.cpp
  lib/packing.cpp
  lib/packer.cpp
  lib/config.cpp
  lib/writers.cpp
  lib/run.cpp)
target_include_directories(randomwalk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(randomwalk PUBLIC poly::common PRIVATE poly_warnings)
add_library(poly::randomwalk ALIAS randomwalk)

add_executable(rw_pack apps/rw_pack.cpp)
target_link_libraries(rw_pack PRIVATE randomwalk poly_warnings)

install(TARGETS rw_pack RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
