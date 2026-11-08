# 3DBFM module: reactive Monte Carlo of radical copolymerization on the 3D bond fluctuation
# model. Builds the static library `bfm` and the command-line front end `bfm_react`.

add_library(bfm STATIC
  lib/bond_set.cpp
  lib/system.cpp
  lib/config.cpp
  lib/writers.cpp
  lib/reaction.cpp
  lib/run.cpp)
target_include_directories(bfm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(bfm PUBLIC poly::common poly::chemistry PRIVATE poly_warnings)
add_library(poly::bfm ALIAS bfm)

add_executable(bfm_react apps/bfm_react.cpp)
target_link_libraries(bfm_react PRIVATE bfm poly_warnings)

install(TARGETS bfm_react RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
