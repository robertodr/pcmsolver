add_library(td_solver-tests OBJECT tdcpcm-iter_gepol-point.cpp 
                                   tdcpcm_gepol-point.cpp 
                                   tdief_gepol-point.cpp 
                                   tdonsagerief_gepol-point.cpp 
                                   tdsingleief_gepol-point.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(td_solver-tests custom_boost)
endif()

# Copy reference files to ${PROJECT_BINARY_DIR}/tests/td_solver (aka ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND reference_files tdcpcm_collocation.npy tdief_lowdin_collocation_K_d.npy
      tdief_lowdin_collocation_tau.npy tdcpcmiterative_collocation.npy   
      tdief_lowdin_collocation_Lambda.npy tdonsagerief_collocation.npy
      tdief_lowdin_collocation_K_0.npy tdief_lowdin_collocation.npy         
      tdsingleief_collocation.npy)
file(COPY ${reference_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# tdcpcm-iter_gepol-point.cpp test
add_Catch_test(tdcpcm-iter_gepol-point "td_solver;tdcpcm;tdcpcm-iter_gepol-point")

# tdcpcm_gepol-point.cpp test
add_Catch_test(tdcpcm_gepol-point "td_solver;tdcpcm;tdcpcm_gepol-point")

# tdsingleief_gepol-point.cpp test
add_Catch_test(tdsingleief_gepol-point "td_solver;tdief;tdsingleief_gepol-point")

# tdonsagerief_gepol-point.cpp test
add_Catch_test(tdonsagerief_gepol-point "td_solver;tdief;tdonsagerief_gepol-point")

# tdief_gepol-point.cpp test
add_Catch_test(tdief_gepol-point "td_solver;tdiefpcm;tdief_gepol-point")

# This executable updates the .npy files containing the reference values
if(BUILD_STANDALONE)
  link_directories(${PROJECT_BINARY_DIR}/lib)
  add_executable(update_tdsolver_reference_files update_tdsolver_reference_files.cpp)
  if(BUILD_CUSTOM_BOOST)
    add_dependencies(update_tdsolver_reference_files custom_boost)
  endif()
  if(STATIC_LIBRARY_ONLY)
    target_link_libraries(update_tdsolver_reference_files pcm-static ${ZLIB_LIBRARIES})
  else()
    target_link_libraries(update_tdsolver_reference_files pcm-shared)
  endif()
endif()
