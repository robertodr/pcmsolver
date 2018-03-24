add_library(td_solver-tests OBJECT
  tdcpcm_gepol-point.cpp
  tdief_gepol-point.cpp
  tdonsagerief_gepol-point.cpp
  tdsingleief_gepol-point.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(td_solver-tests custom_boost)
endif()

# tdcpcm_gepol-point.cpp test
add_Catch_test(
  NAME
    tdcpcm_gepol-point
  LABELS
    td_solvert
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
  target_link_libraries(update_tdsolver_reference_files pcm-static)
  if(BUILD_CUSTOM_BOOST)
    add_dependencies(update_tdsolver_reference_files custom_boost)
  endif()

  add_executable(plot_tdsolver plot_tdsolver.cpp)
  target_link_libraries(plot_tdsolver pcm-static)
  if(BUILD_CUSTOM_BOOST)
    add_dependencies(plot_tdsolver custom_boost)
  endif()
endif()
