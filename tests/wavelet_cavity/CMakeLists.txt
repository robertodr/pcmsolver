set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/wavelet_cavity_NH3.cpp)
add_test(NAME wavelet_cavity_NH3 COMMAND unit_tests.x [wavelet_cavity_NH3])
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/wavelet_cavity_C6H6.cpp)
add_test(NAME wavelet_cavity_C6H6 COMMAND unit_tests.x [wavelet_cavity_C6H6])
