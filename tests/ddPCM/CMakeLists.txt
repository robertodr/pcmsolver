add_library(ddPCM-tests OBJECT ddPCM.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(ddPCM-tests custom_boost)
endif()

# ddPCM.cpp test
add_Catch_test(ddPCM "")

