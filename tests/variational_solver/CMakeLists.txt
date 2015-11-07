# viefpcm_gepol-NH3_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_gepol-NH3_from-file.cpp)
add_Catch_test(viefpcm_gepol-NH3_from-file "variational_solver;viefpcm;viefpcm_gepol-NH3_from-file")
set_tests_properties(viefpcm_gepol-NH3_from-file PROPERTIES DEPENDS viefpcm_gepol-NH3)

# viefpcm_gepol-NH3.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_gepol-NH3.cpp)
add_Catch_test(viefpcm_gepol-NH3 "variational_solver;viefpcm;viefpcm_gepol-NH3")

# viefpcm_anisotropic-gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_anisotropic-gepol-point.cpp)
add_Catch_test(viefpcm_anisotropic-gepol-point "variational_solver;viefpcm;viefpcm_anisotropic-gepol-point;anisotropic")

# viefpcm_diffuse-gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_diffuse-gepol-point.cpp)
add_Catch_test(viefpcm_diffuse-gepol-point "variational_solver;viefpcm;viefpcm_diffuse-gepol-point;anisotropic")

# viefpcm_gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_gepol-point.cpp)
add_Catch_test(viefpcm_gepol-point "variational_solver;viefpcm;viefpcm_gepol-point")

# viefpcm_symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_symmetry.cpp)
add_Catch_test(viefpcm_symmetry "variational_solver;viefpcm;viefpcm_symmetry")

# viefpcm_anisotropic-symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_anisotropic-symmetry.cpp)
add_Catch_test(viefpcm_anisotropic-symmetry "variational_solver;viefpcm;viefpcm_anisotropic-symmetry;anisotropic")

# viefpcm_gepol-C2H4_D2h.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_gepol-C2H4_D2h.cpp)
add_Catch_test(viefpcm_gepol-C2H4_D2h "viefpcm;viefpcm_symmetry;viefpcm_gepol-C2H4_D2h")

# viefpcm_gepol-point_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/viefpcm_gepol-point_from-file.cpp)
add_Catch_test(viefpcm_gepol-point_from-file "variational_solver;viefpcm;viefpcm_gepol-point_from-file")
set_tests_properties(viefpcm_gepol-point_from-file PROPERTIES DEPENDS viefpcm_gepol-point)
