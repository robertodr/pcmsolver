target_sources(pcm-objlib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ITDSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TDCPCMSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TDIEFSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TDSingleIEFSolver.cpp
)

# List of headers
list(APPEND headers_list
  Debye.hpp
  ITDSolver.hpp
  TDCPCMSolver.hpp
  TDIEFSolver.hpp
  TDSingleIEFSolver.hpp
  TDSolver.hpp
  TDSolverData.hpp
)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/td_solver)
endforeach()
