# Code shared by every module: the random number generator and the configuration-file reader.
# Header-only; modules link the interface target poly::common.

add_library(poly_common INTERFACE)
target_include_directories(poly_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(poly::common ALIAS poly_common)

install(DIRECTORY include/ DESTINATION include)
