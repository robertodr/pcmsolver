add_executable(C_host C_host.c)
target_include_directories(C_host
  PRIVATE
    ${PROJECT_SOURCE_DIR}/tests/C_host
  )
target_link_libraries(C_host PCMSolver)
target_link_libraries(C_host m)
set_target_properties(C_host
  PROPERTIES
    LINKER_LANGUAGE C
  )
add_test(
  NAME
    C_host
  COMMAND
    C_host
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(C_host
  PROPERTIES
    LABELS "api;C_host"
  )

add_executable(fail-C_host fail-C_host.c)
target_include_directories(fail-C_host
  PRIVATE
    ${PROJECT_SOURCE_DIR}/tests/C_host
  )
target_link_libraries(fail-C_host PCMSolver)
target_link_libraries(fail-C_host m)
set_target_properties(fail-C_host
  PROPERTIES
    LINKER_LANGUAGE C
  )
add_test(
  NAME
    fail-C_host
  COMMAND
    fail-C_host
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(fail-C_host
  PROPERTIES
    LABELS "api;C_host;fail"
    WILL_FAIL TRUE
  )

add_executable(TD-C_host TD-C_host.c)
target_include_directories(td-C_host
  PRIVATE
    ${PROJECT_SOURCE_DIR}/tests/C_host
    ${PROJECT_SOURCE_DIR}/api
  )
target_link_libraries(TD-C_host PCMSolver)
target_link_libraries(TD-C_host m)
set_target_properties(TD-C_host
  PROPERTIES
    LINKER_LANGUAGE C
  )
file(
  COPY 
    tdsolver.pcm 
  DESTINATION 
    ${CMAKE_CURRENT_BINARY_DIR}
  )
add_test(
  NAME 
    parse-tdsolver 
  COMMAND 
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --inp=tdsolver.pcm
  WORKING_DIRECTORY 
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(parse-tdsolver 
  PROPERTIES 
  LABELS "parser_script"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
add_test(
  NAME
    TD-C_host
  COMMAND
    TD-C_host
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(TD-C_host
  PROPERTIES
    LABELS "api;C_host;TD-C_host;td_solver"
    DEPENDS parse-tdsolver
  )
