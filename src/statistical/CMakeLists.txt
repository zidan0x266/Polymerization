# statistical module: coordinate-free kinetic Monte Carlo (Gillespie) of radical polymerization
# with the shared FRP/ATRP chemistry. Builds the static library `statistical` and the
# command-line front end `stat_react`.

add_library(statistical STATIC
  lib/network.cpp
  lib/kinetics.cpp
  lib/config.cpp
  lib/run.cpp)
target_include_directories(statistical PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(statistical PUBLIC poly::common poly::chemistry PRIVATE poly_warnings)
add_library(poly::statistical ALIAS statistical)

add_executable(stat_react apps/stat_react.cpp)
target_link_libraries(stat_react PRIVATE statistical poly_warnings)

install(TARGETS stat_react RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
