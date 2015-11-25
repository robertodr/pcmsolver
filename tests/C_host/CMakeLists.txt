add_executable(C_host C_host.c)
target_link_libraries(C_host pcm getkw ${ZLIB_LIBRARIES})
set_target_properties(C_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/C_host
                                          INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                          LINKER_LANGUAGE C)
add_test(NAME C_host COMMAND C_host)
set_tests_properties(C_host PROPERTIES LABELS "api;C_host")

add_executable(VPCM-C_host VPCM-C_host.c)
target_link_libraries(VPCM-C_host pcm getkw ${ZLIB_LIBRARIES})
add_custom_command(
  TARGET VPCM-C_host
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/vcpcm.inp ${CMAKE_CURRENT_BINARY_DIR}/pcmsolver.inp
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py ${CMAKE_CURRENT_BINARY_DIR}/pcmsolver.inp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(VPCM-C_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/VPCM-C_host
                                          INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                          LINKER_LANGUAGE C)
add_test(NAME VPCM-C_host COMMAND VPCM-C_host)
set_tests_properties(VPCM-C_host PROPERTIES LABELS "api;C_host;variational_solver")
