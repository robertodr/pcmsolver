if (NOT DEFINED DEFUALT_CXX_FLAGS_SET OR RESET_FLAGS)
  
  if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
      set (CMAKE_CXX_FLAGS "-std=gnu++98 -fPIC") #-std=c++11")
      set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG -Wall -Wextra -Winit-self -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Wwrite-strings ")
      set (CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -DNDEBUG -Wno-unused")
      if (ENABLE_CODE_COVERAGE)
          set (CMAKE_CXX_FLAGS
              "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
          set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
      endif()
  elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
        set (CMAKE_CXX_FLAGS "-std=gnu++98 -fPIC") # -std=c++0x")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -debug -DDEBUG -w2 -wd2196 -Wall -Weffc++ -Wuninitialized")
	set (CMAKE_CXX_FLAGS_RELEASE "-debug -O3 -DNDEBUG -wd2196")
	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
  endif ()
  
  if (CMAKE_CXX_COMPILER_ID MATCHES PGI)
      set(CMAKE_CXX_FLAGS         "-g")
      set(CMAKE_CXX_FLAGS_DEBUG   "-O0")
      set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  endif()
  
  if (CMAKE_CXX_COMPILER_ID MATCHES XL)
      set(CMAKE_CXX_FLAGS         "-g")
      set(CMAKE_CXX_FLAGS_DEBUG   "-O0")
      set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  endif()

  save_compiler_flags(CXX)

endif ()


