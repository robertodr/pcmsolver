configure_file(tdief_gepol-point.cpp.in ${CMAKE_CURRENT_LIST_DIR}/tdief_gepol-point.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/tdief_gepol-point.cpp)
add_test(NAME tdief_gepol-point COMMAND unit_tests.x [tdief_gepol-point])

configure_file(tdcpcm_gepol-point.cpp.in ${CMAKE_CURRENT_LIST_DIR}/tdcpcm_gepol-point.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/tdcpcm_gepol-point.cpp)
add_test(NAME tdcpcm_gepol-point COMMAND unit_tests.x [tdcpcm_gepol-point])

configure_file(tdcpcm-iter_gepol-point.cpp.in ${CMAKE_CURRENT_LIST_DIR}/tdcpcm-iter_gepol-point.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/tdcpcm-iter_gepol-point.cpp)
add_test(NAME tdcpcm-iter_gepol-point COMMAND unit_tests.x [tdcpcm-iter_gepol-point])

configure_file(tdsingleief_gepol-point.cpp.in ${CMAKE_CURRENT_LIST_DIR}/tdsingleief_gepol-point.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/tdsingleief_gepol-point.cpp)
add_test(NAME tdsingleief_gepol-point COMMAND unit_tests.x [tdsingleief_gepol-point])

configure_file(tdonsagerief_gepol-point.cpp.in ${CMAKE_CURRENT_LIST_DIR}/tdonsagerief_gepol-point.cpp @ONLY)
set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/tdonsagerief_gepol-point.cpp)
add_test(NAME tdonsagerief_gepol-point COMMAND unit_tests.x [tdonsagerief_gepol-point])

list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})
add_executable(update_tdsolver_reference_files.x update_tdsolver_reference_files.cpp)
target_link_libraries(update_tdsolver_reference_files.x pcm ${external_libraries})
