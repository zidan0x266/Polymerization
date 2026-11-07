# dll module: the dynamic lattice liquid model on a fully occupied FCC lattice, with the shared
# FRP/ATRP chemistry. Builds the static library `dll` and the command-line front end `dll_react`.

add_library(dll STATIC
  lib/system.cpp
  lib/motion.cpp
  lib/crossing.cpp
  lib/reaction.cpp
  lib/initial.cpp
  lib/config.cpp
  lib/writers.cpp
  lib/run.cpp)
target_include_directories(dll PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(dll PUBLIC poly::common poly::chemistry PRIVATE poly::randomwalk poly_warnings)
add_library(poly::dll ALIAS dll)

add_executable(dll_react apps/dll_react.cpp)
target_link_libraries(dll_react PRIVATE dll poly_warnings)

install(TARGETS dll_react RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
