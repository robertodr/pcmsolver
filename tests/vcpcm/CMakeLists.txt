# vcpcm_gepol-point_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm_gepol-point_from-file.cpp)
add_Catch_test(vcpcm_gepol-point_from-file "variational_solver;vcpcm;vcpcm_gepol-point_from-file")
set_tests_properties(vcpcm_gepol-point_from-file PROPERTIES DEPENDS vcpcm_gepol-point)

# vcpcm_gepol-NH3_from-file.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm_gepol-NH3_from-file.cpp)
add_Catch_test(vcpcm_gepol-NH3_from-file "variational_solver;vcpcm;vcpcm_gepol-NH3_from-file")
set_tests_properties(vcpcm_gepol-NH3_from-file PROPERTIES DEPENDS vcpcm_gepol-NH3)

# vcpcm-guess_gepol-C2H4_D2h.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm-guess_gepol-C2H4_D2h.cpp)
add_Catch_test(vcpcm-guess_gepol-C2H4_D2h "variational_solver;initial_guess;vcpcm;vcpcm-guess_gepol-C2H4_D2h")

# vcpcm_gepol-C2H4_D2h.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm_gepol-C2H4_D2h.cpp)
add_Catch_test(vcpcm_gepol-C2H4_D2h "variational_solver;vcpcm;vcpcm_symmetry;vcpcm_gepol-C2H4_D2h")

# vcpcm_gepol-point.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm_gepol-point.cpp)
add_Catch_test(vcpcm_gepol-point "variational_solver;vcpcm;vcpcm_gepol-point")

# vcpcm_symmetry.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm_symmetry.cpp)
add_Catch_test(vcpcm_symmetry "variational_solver;vcpcm;vcpcm_symmetry")

# vcpcm_gepol-NH3.cpp test
set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/vcpcm_gepol-NH3.cpp)
add_Catch_test(vcpcm_gepol-NH3 "variational_solver;vcpcm;vcpcm_gepol-NH3")
