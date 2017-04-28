/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
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

#include "catch.hpp"

#include <iostream>

#include <Eigen/Core>

#include "utils/Molecule.hpp"
#include "TestingMolecules.hpp"
#include "solver/ddPCM.hpp"

using namespace pcm;
using solver::ddPCM;

/*! \class ddPCM
 *
 */
TEST_CASE("ddCOSMO solver with NH3 molecule", "[ddPCM]") {
  Molecule molec = NH3();
  ddPCM solver(molec);
  Eigen::VectorXd potential = computeMEP(molec, solver.cavity());
  Eigen::MatrixXd charges = solver.computeCharges(potential);
}
