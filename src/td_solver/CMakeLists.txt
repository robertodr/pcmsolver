# List of headers
list(APPEND headers_list Debye.hpp
                         TDSolver.hpp
                         TDCPCMSolver.hpp
                         TDIEFSolver.hpp
                         ITDSolver.hpp
                         TDSingleIEFSolver.hpp
                         TDSolverData.hpp)

# List of sources
list(APPEND sources_list ITDSolver.cpp
                         TDIEFSolver.cpp
                         TDCPCMSolver.cpp
                         TDSingleIEFSolver.cpp)

add_library(td_solver OBJECT ${sources_list} ${headers_list})
set_target_properties(td_solver PROPERTIES POSITION_INDEPENDENT_CODE 1 )
add_dependencies(td_solver generate-config-hpp)
target_compile_options(td_solver PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/td_solver)
endforeach()
