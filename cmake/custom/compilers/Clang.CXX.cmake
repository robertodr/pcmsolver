# Overrides contents of all variables previously set by CMake
if(NOT DEFINED ENV{CXXFLAGS})
    if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
        if(HAS_CXX11_SUPPORT)
            set(CMAKE_CXX_FLAGS "-fPIC ${CXX11_COMPILER_FLAGS}")
        else()
            set(CMAKE_CXX_FLAGS "-fPIC -std=gnu++98")
        endif()
        set(CMAKE_CXX_FLAGS_DEBUG    "-O0 -DDEBUG -Wall -Wextra -Winit-self -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Wwrite-strings -Weffc++ -Wdocumentation")
        set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG -Wno-unused")
    endif()
endif()
