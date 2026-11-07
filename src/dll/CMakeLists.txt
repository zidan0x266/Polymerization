# dll module: the dynamic lattice liquid model on a fully occupied FCC lattice. Builds the static
# library `dll` and the command-line front end `dll_react`. Phase 1: motion of chains and solvent.

add_library(dll STATIC
  lib/system.cpp
  lib/motion.cpp
  lib/initial.cpp
  lib/config.cpp
  lib/writers.cpp
  lib/run.cpp)
target_include_directories(dll PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(dll PUBLIC poly::common PRIVATE poly::randomwalk poly_warnings)
add_library(poly::dll ALIAS dll)

add_executable(dll_react apps/dll_react.cpp)
target_link_libraries(dll_react PRIVATE dll poly_warnings)

install(TARGETS dll_react RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
