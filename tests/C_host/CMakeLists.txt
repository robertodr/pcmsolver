link_directories(${PROJECT_BINARY_DIR}/lib)
add_executable(C_host C_host.c)
add_executable(fail-C_host fail-C_host.c)
add_executable(TD-C_host TD-C_host.c)
if(STATIC_LIBRARY_ONLY)
  target_link_libraries(C_host pcm-static ${ZLIB_LIBRARIES})
  target_link_libraries(fail-C_host pcm-static ${ZLIB_LIBRARIES})
  target_link_libraries(TD-C_host pcm-static ${ZLIB_LIBRARIES})
else()
  target_link_libraries(C_host pcm-shared)
  target_link_libraries(fail-C_host pcm-shared)
  target_link_libraries(TD-C_host pcm-shared)
endif()

target_link_libraries(C_host m)
set_target_properties(C_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/C_host
                                        INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                        LINKER_LANGUAGE C)
add_test(NAME C_host COMMAND C_host
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(C_host PROPERTIES LABELS "api;C_host")

target_link_libraries(fail-C_host m)
set_target_properties(fail-C_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/C_host
                                        INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                        LINKER_LANGUAGE C)
add_test(NAME fail-C_host COMMAND fail-C_host
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(fail-C_host PROPERTIES LABELS "api;C_host;fail")
set_tests_properties(fail-C_host PROPERTIES WILL_FAIL TRUE)

target_link_libraries(TD-C_host m)
set_target_properties(TD-C_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/C_host
                                        INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                        LINKER_LANGUAGE C)
file(COPY tdsolver.pcm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/tdsolver.pcm ${CMAKE_CURRENT_BINARY_DIR}/pcmsolver.inp)
add_test(NAME parse-tdsolver COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/pcmsolver.py pcmsolver.inp
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(parse-tdsolver PROPERTIES LABELS "parser_script" DEPENDS generate-pcmsolver-py)
add_test(NAME TD-C_host
         COMMAND TD-C_host
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(TD-C_host PROPERTIES LABELS "api;C_host;tdsolver"
                                          DEPENDS parse-tdsolver)
