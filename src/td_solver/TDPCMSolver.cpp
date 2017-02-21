/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#include "TDPCMSolver.hpp"

#include <string>
#include <sstream>

#include "Config.hpp"

#include "Debye.hpp"

TDPCMSolver::TDPCMSolver(double es, double ed, double t)
    : permittivity_(Debye(es, ed, t)), built_(false) {}

std::string TDPCMSolver::printEnvironment() {
  std::stringstream tmp;
  tmp << ".... Inside " << std::endl;
  tmp << "Green's function type: vacuum" << std::endl;
  tmp << ".... Outside " << std::endl;
  tmp << "Green's function type: uniform dielectric" << std::endl;
  tmp << permittivity_;
  return tmp.str();
}
