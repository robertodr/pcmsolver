#!/usr/bin/env python

# This file is autogenerated by Autocmake v1.0.0-alpha-x http://autocmake.org
# Copyright (c) 2015-2016 by Radovan Bast, Jonas Juselius, and contributors.

import os
import sys

sys.path.insert(0, 'cmake')
from autocmake import configure
from autocmake.external import docopt


options = """
Usage:
  ./setup.py [options] [<builddir>]
  ./setup.py (-h | --help)

Options:
  --fc=<FC>                                Fortran compiler [default: gfortran].
  --extra-fc-flags=<EXTRA_FCFLAGS>         Extra Fortran compiler flags [default: ''].
  --cc=<CC>                                C compiler [default: gcc].
  --extra-cc-flags=<EXTRA_CFLAGS>          Extra C compiler flags [default: ''].
  --cxx=<CXX>                              C++ compiler [default: g++].
  --extra-cxx-flags=<EXTRA_CXXFLAGS>       Extra C++ compiler flags [default: ''].
  --bindir=<CMAKE_INSTALL_BINDIR>          User executables [default: bin].
  --libdir=<CMAKE_INSTALL_LIBDIR>          Object code libraries [default: lib].
  --includedir=<CMAKE_INSTALL_INCLUDEDIR>  C header files [default: include].
  --datadir=<CMAKE_INSTALL_DATADIR>        Read-only architecture-independent data root [default: share].
  --ccache=<USE_CCACHE>                    Toggle use of ccache <ON/OFF> [default: ON].
  --add-definitions=<STRING>               Add preprocesor definitions [default: ''].
  --coverage                               Enable code coverage [default: False].
  --int64                                  Enable 64bit integers [default: False].
  --omp                                    Enable OpenMP parallelization [default: False].
  --python=<PYTHON_INTERPRETER>            The Python interpreter (development version) to use. [default: ''].
  --fbindings=<ENABLE_FORTRAN_API>         Enable compilation of Fortran 90 API bindings <ON/OFF> [default: ON].
  --boost-headers=<BOOST_INCLUDEDIR>       Include directories for Boost [default: ''].
  --boost-libraries=<BOOST_LIBRARYDIR>     Library directories for Boost [default: ''].
  --build-boost=<FORCE_CUSTOM_BOOST>       Deactivate Boost detection and build on-the-fly <ON/OFF> [default: OFF].
  --static                                 Create only the static library [default: False].
  --eigen=<EIGEN3_ROOT>                    Root directory for Eigen3 [default: ''].
  --type=<TYPE>                            Set the CMake build type (debug, release, or relwithdeb) [default: release].
  --generator=<STRING>                     Set the CMake build system generator [default: Unix Makefiles].
  --show                                   Show CMake command and exit.
  --cmake-executable=<CMAKE_EXECUTABLE>    Set the CMake executable [default: cmake].
  --cmake-options=<STRING>                 Define options to CMake [default: ''].
  --prefix=<PATH>                          Set the install path for make install.
  <builddir>                               Build directory.
  -h --help                                Show this screen.
"""


def gen_cmake_command(options, arguments):
    """
    Generate CMake command based on options and arguments.
    """
    command = []
    command.append('FC={0}'.format(arguments['--fc']))
    command.append('CC={0}'.format(arguments['--cc']))
    command.append('CXX={0}'.format(arguments['--cxx']))
    command.append(arguments['--cmake-executable'])
    command.append('-DEXTRA_FCFLAGS="{0}"'.format(arguments['--extra-fc-flags']))
    command.append('-DEXTRA_CFLAGS="{0}"'.format(arguments['--extra-cc-flags']))
    command.append('-DEXTRA_CXXFLAGS="{0}"'.format(arguments['--extra-cxx-flags']))
    command.append('-DCMAKE_INSTALL_BINDIR={0}'.format(arguments['--bindir']))
    command.append('-DCMAKE_INSTALL_LIBDIR={0}'.format(arguments['--libdir']))
    command.append('-DCMAKE_INSTALL_INCLUDEDIR={0}'.format(arguments['--includedir']))
    command.append('-DCMAKE_INSTALL_DATADIR={0}'.format(arguments['--datadir']))
    command.append('-DUSE_CCACHE={0}'.format(arguments['--ccache']))
    command.append('-DPREPROCESSOR_DEFINITIONS="{0}"'.format(arguments['--add-definitions']))
    command.append('-DENABLE_CODE_COVERAGE={0}'.format(arguments['--coverage']))
    command.append('-DENABLE_64BIT_INTEGERS={0}'.format(arguments['--int64']))
    command.append('-DENABLE_OPENMP={0}'.format(arguments['--omp']))
    command.append('-DPYTHON_INTERPRETER="{0}"'.format(arguments['--python']))
    command.append('-DENABLE_FORTRAN_API={0}'.format(arguments['--fbindings']))
    command.append('-DBOOST_INCLUDEDIR="{0}"'.format(arguments['--boost-headers']))
    command.append('-DBOOST_LIBRARYDIR="{0}"'.format(arguments['--boost-libraries']))
    command.append('-DFORCE_CUSTOM_BOOST={0}'.format(arguments['--build-boost']))
    command.append('-DBOOST_MINIMUM_REQUIRED="1.54.0"')
    command.append('-DBOOST_COMPONENTS_REQUIRED=""')
    command.append('-DSTATIC_LIBRARY_ONLY={0}'.format(arguments['--static']))
    command.append('-DEIGEN3_ROOT="{0}"'.format(arguments['--eigen']))
    command.append('-DCMAKE_BUILD_TYPE={0}'.format(arguments['--type']))
    command.append('-G "{0}"'.format(arguments['--generator']))
    if arguments['--cmake-options'] != "''":
        command.append(arguments['--cmake-options'])
    if arguments['--prefix']:
        command.append('-DCMAKE_INSTALL_PREFIX="{0}"'.format(arguments['--prefix']))

    return ' '.join(command)


# parse command line args
try:
    arguments = docopt.docopt(options, argv=None)
except docopt.DocoptExit:
    sys.stderr.write('ERROR: bad input to {0}\n'.format(sys.argv[0]))
    sys.stderr.write(options)
    sys.exit(-1)


# use extensions to validate/post-process args
if configure.module_exists('extensions'):
    import extensions
    arguments = extensions.postprocess_args(sys.argv, arguments)


root_directory = os.path.dirname(os.path.realpath(__file__))


build_path = arguments['<builddir>']


# create cmake command
cmake_command = '{0} {1}'.format(gen_cmake_command(options, arguments), root_directory)


# run cmake
configure.configure(root_directory, build_path, cmake_command, arguments['--show'])
