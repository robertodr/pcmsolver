# Executables for testing separate from Boost unit test framework or for profiling purposes
# The executables will be created in ${PROJECT_BINARY_DIR}/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

list(APPEND external_libraries ${ZLIB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY})

#add_executable(debug_wavcav.x debug_wavcav.cpp)
#target_link_libraries(debug_wavcav.x solver cavity wavcav pwl wem utils)

#add_executable(plot_green_spherical.x plot_green_spherical.cpp)
#target_link_libraries(plot_green_spherical.x utils ${external_libraries})

#add_executable(plot_green_spherical-CASE1.x plot_green_spherical-CASE1.cpp)
#target_link_libraries(plot_green_spherical-CASE1.x utils ${external_libraries})
#add_executable(plot_green_spherical-CASE2.x plot_green_spherical-CASE2.cpp)
#target_link_libraries(plot_green_spherical-CASE2.x utils ${external_libraries})
#add_executable(plot_green_spherical-CASE3.x plot_green_spherical-CASE3.cpp)
#target_link_libraries(plot_green_spherical-CASE3.x utils ${external_libraries})
#add_executable(plot_green_spherical-CASE4.x plot_green_spherical-CASE4.cpp)
#target_link_libraries(plot_green_spherical-CASE4.x utils ${external_libraries})

#add_executable(check_Coulomb_coefficient.x check_Coulomb_coefficient.cpp)
#target_link_libraries(check_Coulomb_coefficient.x utils ${external_libraries})

add_executable(run_pcm.x run_pcm.cpp)
if(BUILD_CUSTOM_BOOST)
    add_dependencies(run_pcm.x custom_boost)
endif()
target_link_libraries(run_pcm.x pcm getkw ${external_libraries})
set_property(GLOBAL PROPERTY PCMSolver_EXECUTABLE ${PROJECT_BINARY_DIR}/bin/run_pcm.x)

add_executable(profile_pwc.x profile_pwc.cpp)
target_link_libraries(profile_pwc.x pcm getkw ${external_libraries})

add_executable(profile_pwl.x profile_pwl.cpp)
target_link_libraries(profile_pwl.x pcm getkw ${external_libraries})

add_executable(benzene_wavelet.x benzene_wavelet.cpp)
target_link_libraries(benzene_wavelet.x pcm getkw ${external_libraries})
#add_reference(${CMAKE_SOURCE_DIR}/examples/benzene2.dat ${CMAKE_CURRENT_BINARY_DIR})
#add_reference(${CMAKE_SOURCE_DIR}/examples/sphere2.dat ${CMAKE_CURRENT_BINARY_DIR})

add_executable(poly_wavelet.x poly_wavelet.cpp)
target_link_libraries(poly_wavelet.x pcm getkw ${external_libraries})
#add_reference(${CMAKE_SOURCE_DIR}/examples/C32H66_lev2.dat ${CMAKE_CURRENT_BINARY_DIR})
#add_reference(${CMAKE_SOURCE_DIR}/examples/C8H18_lev2.dat ${CMAKE_CURRENT_BINARY_DIR})

add_executable(poly_linear.x poly_linear.cpp)
target_link_libraries(poly_linear.x pcm getkw ${external_libraries})

add_executable(molecule_wavcav_wavelet.x molecule_wavcav_wavelet.cpp)
target_link_libraries(molecule_wavcav_wavelet.x pcm getkw ${external_libraries})

unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

