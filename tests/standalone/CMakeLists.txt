# Copy input files to ${PROJECT_BINARY_DIR}/tests/standalone (aka ${CMAKE_CURRENT_BINARY_DIR})
file(
  COPY
    onsager_tdsolver.pcm
    standalone.pcm
    standalone_atoms.pcm
    standalone_explicit.pcm
    standalone_tdsolver.pcm
  DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}
  )

add_test(
  NAME
    pcmsolver_standalone
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone
  PROPERTIES
  LABELS "pcmsolver;standalone"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-atoms
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_atoms.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-atoms
  PROPERTIES
  LABELS "pcmsolver;standalone;atoms"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-explicit
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_explicit.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-explicit
  PROPERTIES
  LABELS "pcmsolver;standalone;explicit"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-tdsolver
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=standalone_tdsolver.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-tdsolver
  PROPERTIES
  LABELS "pcmsolver;standalone;explicit"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )

add_test(
  NAME
    pcmsolver_standalone-onsager_tdsolver
  COMMAND
    ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/go_pcm.py
      --exe=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
      --inp=onsager_tdsolver.pcm
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(pcmsolver_standalone-onsager_tdsolver
  PROPERTIES
  LABELS "pcmsolver;standalone;onsager_tdsolver"
  ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
  )
