#!/usr/bin/env python

# This file is autogenerated by Autocmake http://autocmake.org
# Copyright (c) 2015 by Radovan Bast and Jonas Juselius

import os
import sys

sys.path.append('cmake/lib/docopt')

sys.path.append('cmake/lib')
from config import configure
import docopt


options = """
Usage:
  ./setup.py [options] [<builddir>]
  ./setup.py (-h | --help)

Options:
  --cxx=<CXX>                            C++ compiler [default: g++].
  --extra-cxx-flags=<EXTRA_CXXFLAGS>     Extra C++ compiler flags [default: ''].
  --cc=<CC>                              C compiler [default: gcc].
  --extra-cc-flags=<EXTRA_CFLAGS>        Extra C compiler flags [default: ''].
  --fc=<FC>                              Fortran compiler [default: gfortran].
  --extra-fc-flags=<EXTRA_FCFLAGS>       Extra Fortran compiler flags [default: ''].
  --coverage                             Enable code coverage [default: False].
  --int64                                Enable 64bit integers [default: False].
  --omp                                  Enable OpenMP parallelization [default: False].
  --python=<PYTHON_INTERPRETER>          The Python interpreter (development version) to use. [default: ''].
  --boost-headers=<BOOST_INCLUDEDIR>     Include directories for Boost [default: ''].
  --boost-libraries=<BOOST_LIBRARYDIR>   Library directories for Boost [default: ''].
  --build-boost=<FORCE_CUSTOM_BOOST>     Deactivate Boost detection and build on-the-fly <ON/OFF> [default: OFF].
  --type=<TYPE>                          Set the CMake build type (debug, release, or relwithdeb) [default: release].
  --generator=<STRING>                   Set the CMake build system generator [default: Unix Makefiles].
  --show                                 Show CMake command and exit.
  --cmake-executable=<CMAKE_EXECUTABLE>  Set the CMake executable [default: cmake].
  --cmake-options=<STRING>               Define options to CMake [default: ''].
  <builddir>                             Build directory.
  -h --help                              Show this screen.
"""


def gen_cmake_command(options, arguments):
    """
    Generate CMake command based on options and arguments.
    """
    command = []
    command.append('CXX=%s' % arguments['--cxx'])
    command.append('CC=%s' % arguments['--cc'])
    command.append('FC=%s' % arguments['--fc'])
    command.append('%s' % arguments['--cmake-executable'])
    command.append('-DEXTRA_CXXFLAGS="%s"' % arguments['--extra-cxx-flags'])
    command.append('-DEXTRA_CFLAGS="%s"' % arguments['--extra-cc-flags'])
    command.append('-DEXTRA_FCFLAGS="%s"' % arguments['--extra-fc-flags'])
    command.append('-DENABLE_CODE_COVERAGE=%s' % arguments['--coverage'])
    command.append('-DENABLE_64BIT_INTEGERS=%s' % arguments['--int64'])
    command.append('-DENABLE_OPENMP=%s' % arguments['--omp'])
    command.append('-DPYTHON_INTERPRETER="%s"' % arguments['--python'])
    command.append('-DBOOST_INCLUDEDIR="%s"' % arguments['--boost-headers'])
    command.append('-DBOOST_LIBRARYDIR="%s"' % arguments['--boost-libraries'])
    command.append('-DFORCE_CUSTOM_BOOST="%s"' % arguments['--build-boost'])
    command.append('-DCMAKE_BUILD_TYPE=%s' % arguments['--type'])
    command.append('-G "%s"' % arguments['--generator'])
    if arguments['--cmake-options'] != "''":
        command.append('%s' % arguments['--cmake-options'])

    return ' '.join(command)


try:
    arguments = docopt.docopt(options, argv=None)
except docopt.DocoptExit:
    sys.stderr.write('ERROR: bad input to %s\n' % sys.argv[0])
    sys.stderr.write(options)
    sys.exit(-1)

root_directory = os.path.dirname(os.path.realpath(__file__))
build_path = arguments['<builddir>']
cmake_command = '%s %s' % (gen_cmake_command(options, arguments), root_directory)
configure(root_directory, build_path, cmake_command, arguments['--show'])
