# Take care of updating the cache for fresh configurations
if (NOT DEFINED HAVE_CXX_FLAGS)
  if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-O2 -Wall -DNDEBUG -Wno-unused")
	if (ENABLE_CODE_COVERAGE)
	  set (CMAKE_CXX_FLAGS 
		"${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
	  set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
	endif()
  elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
	set (CMAKE_CXX_FLAGS "-w -Wno-unknown-pragmas")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -debug -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-debug -O3 -DNDEBUG")
	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
  endif ()
  SaveCompilerFlags(CXX)
endif ()


