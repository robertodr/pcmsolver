# List of headers
list(APPEND headers_list VCPCMSolver.hpp VIEFSolver.hpp VPCMSolver.hpp VSolverImpl.hpp RegisterVSolverToFactory.hpp)

# List of sources
list(APPEND sources_list VPCMSolver.cpp VCPCMSolver.cpp VIEFSolver.cpp)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION include/variational_solver)
endforeach()
