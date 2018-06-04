target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tdcpcm_gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tdief_gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tdonsagerief_gepol-point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tdsingleief_gepol-point.cpp
  )

# tdcpcm_gepol-point.cpp test
add_Catch_test(
  NAME
    tdcpcm_gepol-point
  LABELS
    td_solver
    tdcpcm
    tdcpcm_gepol-point
  REFERENCE_FILES
    tdcpcm_collocation.npy
  )

# tdsingleief_gepol-point.cpp test
add_Catch_test(
  NAME
    tdsingleief_gepol-point
  LABELS
    td_solver
    tdief
    tdsingleief_gepol-point
  REFERENCE_FILES
    tdsingleief_collocation.npy
  )

# tdonsagerief_gepol-point.cpp test
add_Catch_test(
  NAME
    tdonsagerief_gepol-point
  LABELS
    td_solver
    tdief
    tdonsagerief_gepol-point
  REFERENCE_FILES
    tdonsagerief_collocation.npy
  )

# tdief_gepol-point.cpp test
add_Catch_test(
  NAME
    tdief_gepol-point
  LABELS
    td_solver
    tdiefpcm
    tdief_gepol-point
  REFERENCE_FILES
    tdief_lowdin_collocation.npy
    tdief_lowdin_collocation_K_0.npy
    tdief_lowdin_collocation_K_d.npy
    tdief_lowdin_collocation_Lambda.npy
    tdief_lowdin_collocation_tau.npy
  )

# This executable updates the .npy files containing the reference values
if(BUILD_STANDALONE)
  add_executable(update_tdsolver_reference_files update_tdsolver_reference_files.cpp)
  target_compile_definitions(update_tdsolver_reference_files
    PRIVATE
      "$<$<BOOL:${HAS_CXX11}>:HAS_CXX11>"
    )
  target_link_libraries(update_tdsolver_reference_files pcm-static)
  target_include_directories(update_tdsolver_reference_files
    PRIVATE
      ${PROJECT_SOURCE_DIR}/tests
    )

  add_executable(plot_tdsolver plot_tdsolver.cpp)
  target_compile_definitions(plot_tdsolver
    PRIVATE
      "$<$<BOOL:${HAS_CXX11}>:HAS_CXX11>"
    )
  target_link_libraries(plot_tdsolver pcm-static)
endif()
