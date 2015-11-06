#!/usr/bin/env python

# This file is autogenerated by Autocmake http://autocmake.org
# Copyright (c) 2015 by Radovan Bast and Jonas Juselius

import os
import sys

sys.path.insert(0, 'cmake')
sys.path.insert(0, 'cmake/lib')
sys.path.insert(0, 'cmake/lib/docopt')
import config
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
  --exdiag                               Enable C++ extended diagnostics flags [default: False].
  --ccache=<USE_CCACHE>                  Toggle use of ccache <ON/OFF> [default: ON].
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
    command.append('CXX={0}'.format(arguments['--cxx']))
    command.append('CC={0}'.format(arguments['--cc']))
    command.append('FC={0}'.format(arguments['--fc']))
    command.append('%s' % arguments['--cmake-executable'])
    command.append('-DEXTRA_CXXFLAGS="{0}"'.format(arguments['--extra-cxx-flags']))
    command.append('-DEXTRA_CFLAGS="{0}"'.format(arguments['--extra-cc-flags']))
    command.append('-DEXTRA_FCFLAGS="{0}"'.format(arguments['--extra-fc-flags']))
    command.append('-DENABLE_EXTENDED_DIAGNOSTICS=%s' % arguments['--exdiag'])
    command.append('-DUSE_CCACHE="{0}"'.format(arguments['--ccache']))
    command.append('-DENABLE_CODE_COVERAGE=%s' % arguments['--coverage'])
    command.append('-DENABLE_64BIT_INTEGERS=%s' % arguments['--int64'])
    command.append('-DENABLE_OPENMP=%s' % arguments['--omp'])
    command.append('-DPYTHON_INTERPRETER="%s"' % arguments['--python'])
    command.append('-DBOOST_INCLUDEDIR="{0}"'.format(arguments['--boost-headers']))
    command.append('-DBOOST_LIBRARYDIR="{0}"'.format(arguments['--boost-libraries']))
    command.append('-DFORCE_CUSTOM_BOOST="{0}"'.format(arguments['--build-boost']))
    command.append('-DBOOST_MINIMUM_REQUIRED="1.54.0"')
    command.append('-DBOOST_COMPONENTS_REQUIRED="''"')
    command.append('-DCMAKE_BUILD_TYPE=%s' % arguments['--type'])
    command.append('-G "%s"' % arguments['--generator'])
    if arguments['--cmake-options'] != "''":
        command.append('%s' % arguments['--cmake-options'])

    return ' '.join(command)


# parse command line args
try:
    arguments = docopt.docopt(options, argv=None)
except docopt.DocoptExit:
    sys.stderr.write('ERROR: bad input to %s\n' % sys.argv[0])
    sys.stderr.write(options)
    sys.exit(-1)


# use extensions to validate/post-process args
if config.module_exists('extensions'):
    import extensions
    arguments = extensions.postprocess_args(sys.argv, arguments)


root_directory = os.path.dirname(os.path.realpath(__file__))


build_path = arguments['<builddir>']


# create cmake command
cmake_command = '%s %s' % (gen_cmake_command(options, arguments), root_directory)


# run cmake
config.configure(root_directory, build_path, cmake_command, arguments['--show'])
