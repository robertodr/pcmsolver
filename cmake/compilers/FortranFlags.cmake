if (NOT DEFINED DEFAULT_Fortran_FLAGS_SET OR RESET_FLAGS)

  if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
  	set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -DVAR_GFORTRAN -DGFORTRAN=445 -fimplicit-none -fPIC -fautomatic -std=f2003")
  	set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -fbacktrace -Wall -Wextra -fcheck=all")
  	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -funroll-all-loops -ftree-vectorize")
        if(ENABLE_64BIT_INTEGERS)                                                 
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
        endif()
        if(ENABLE_BOUNDS_CHECK)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check")
        endif()
        if(ENABLE_CODE_COVERAGE)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
        endif()
        if(ENABLE_VECTORIZATION)                                                    
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_ARCHITECTURE_FLAGS}")
        endif()	      
  endif()
  
  if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
      add_definitions(-DVAR_IFORT)
      set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -w -fpp -assume byterecl -traceback -fPIC -nosave")
      set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -warn all -ftrapuv -check all")
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -ip")
      if(ENABLE_64BIT_INTEGERS)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
      endif()
      if(ENABLE_BOUNDS_CHECK)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion")
      endif()
      if(ENABLE_VECTORIZATION)                                                    
      	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_ARCHITECTURE_FLAGS}")
      endif()	      
  endif()
  
  if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
  	set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wno=155 -fno-second-underscore -DVAR_G95 -fsloppy-char -fPIC")
  	set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -ftrace=full")
  	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
      if(ENABLE_64BIT_INTEGERS)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
      endif()
      if(ENABLE_BOUNDS_CHECK)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check")
      endif()
      if(ENABLE_CODE_COVERAGE)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
      endif()
  endif()
  
  
  if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  	set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -DVAR_PGF90")
  	set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -g -O0 -Mframe")
  	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2 -mcmodel=medium -fast -Munroll")
      if(ENABLE_64BIT_INTEGERS)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -i8")
      else()
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32")
      endif()
      if(ENABLE_BOUNDS_CHECK)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
      endif()
      if(ENABLE_CODE_COVERAGE)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
      endif()
  endif()
  
  if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
  	set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -qzerosize -qextname")
  	set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -g")
  	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2")
      if(ENABLE_64BIT_INTEGERS)
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -q64")
      endif()
  endif()
  
  save_compiler_flags(Fortran)


endif()
