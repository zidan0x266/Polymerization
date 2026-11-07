# Tests for every module. Unit tests are plain executables returning non-zero on failure;
# smoke tests run the real binaries on the small systems in examples/.

set(POLY_EXAMPLES ${CMAKE_SOURCE_DIR}/examples)

# ---- 3dbfm --------------------------------------------------------------------------------------
add_executable(bfm_tests 3dbfm/test_main.cpp)
target_link_libraries(bfm_tests PRIVATE poly::bfm poly_warnings)
add_test(NAME 3dbfm.unit COMMAND bfm_tests)

add_test(NAME 3dbfm.smoke_atrp
  COMMAND bfm_react ${POLY_EXAMPLES}/example_1/smoke_small.ini mechanism=atrp k_p1=0.5 k_p2=0.25
          directory=${CMAKE_BINARY_DIR}/smoke/atrp seed=7 stop_conversion=0.30 snapshot_start=0.10 snapshot_step=0.10)
add_test(NAME 3dbfm.smoke_frp
  COMMAND bfm_react ${POLY_EXAMPLES}/example_1/smoke_small.ini mechanism=frp k_init=0.01 k_p1=0.5 k_p2=0.25 k_term=0.95
          directory=${CMAKE_BINARY_DIR}/smoke/frp seed=11 stop_conversion=0.30 snapshot_start=0.10 snapshot_step=0.10)

# ---- randomwalk ---------------------------------------------------------------------------------
add_executable(rw_tests randomwalk/test_main.cpp)
target_link_libraries(rw_tests PRIVATE poly::randomwalk poly_warnings)
add_test(NAME randomwalk.unit COMMAND rw_tests)

add_test(NAME randomwalk.smoke
  COMMAND rw_pack ${POLY_EXAMPLES}/example_2/smoke_small.ini directory=${CMAKE_BINARY_DIR}/smoke/randomwalk seed=3)

# ---- dll ----------------------------------------------------------------------------------------
add_executable(dll_tests dll/test_main.cpp)
target_link_libraries(dll_tests PRIVATE poly::dll poly_warnings)
add_test(NAME dll.unit COMMAND dll_tests)

add_test(NAME dll.smoke
  COMMAND dll_react ${POLY_EXAMPLES}/example_3/smoke_small.ini directory=${CMAKE_BINARY_DIR}/smoke/dll seed=3)
add_test(NAME dll.smoke_atrp
  COMMAND dll_react ${POLY_EXAMPLES}/example_3/smoke_react.ini mechanism=atrp directory=${CMAKE_BINARY_DIR}/smoke/dll_atrp seed=5)
add_test(NAME dll.smoke_frp
  COMMAND dll_react ${POLY_EXAMPLES}/example_3/smoke_react.ini mechanism=frp reaction_partner=all_neighbors
          directory=${CMAKE_BINARY_DIR}/smoke/dll_frp seed=7)
