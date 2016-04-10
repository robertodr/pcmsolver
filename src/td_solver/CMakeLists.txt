# List of headers
list(APPEND headers_list Debye.hpp RegisterTDSolverToFactory.hpp
                         TDCPCMIterativeSolver.hpp TDCPCMSolver.hpp
                         TDIEFSolver.hpp TDOnsagerIEFSolver.hpp
                         TDPCMIterativeSolver.hpp TDPCMSolver.hpp
                         TDSingleIEFSolver.hpp TDSolverHelperFunctions.hpp
                         TDSolverData.hpp)

# List of sources
list(APPEND sources_list )

add_library(td_solver OBJECT ${sources_list} ${headers_list})
set_target_properties(td_solver PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(td_solver generate-config-hpp)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION include/td_solver)
endforeach()
