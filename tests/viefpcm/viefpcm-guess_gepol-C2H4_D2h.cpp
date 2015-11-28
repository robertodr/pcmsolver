/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include <iostream>

#include "Config.hpp"

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "Vacuum.hpp"
#include "TestingMolecules.hpp"
#include "UniformDielectric.hpp"
#include "VIEFSolver.hpp"
#include "Symmetry.hpp"

/*! \class VIEFSolver
 *  \test \b C2H4D2hGePol tests initial ASC guesses for VIEFSolver using C2H4 in D2h and a GePol cavity
 */
TEST_CASE("Test variational solver initial ASC guesses for the IEFPCM with C2H4 molecule and a GePol cavity",
          "[variational_solver][initial_guess][viefpcm][viefpcm-guess_gepol-C2H4_D2h]")
{
  Molecule molec = C2H4_D2h();
  double area = 10.0;
  double probeRadius = 10.0;
  double minRadius = 100.0;
  GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius);

  double permittivity = 78.39;
  Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
  UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
    UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);

  double Ccharge = 6.0;
  double Hcharge = 1.0;
  double total_charge = 2.0 * Ccharge + 4.0 * Hcharge;
  size_t size = cavity.size();
  int nr_irrep = cavity.pointGroup().nrIrrep();
  Eigen::VectorXd fake_mep = computeMEP(molec, cavity.elements());

  SECTION("Trivial initial guess")
  {
    VIEFSolver solver(VPCMSolver::Trivial, VPCMSolver::SSD);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);
    Eigen::VectorXd guess = solver.initialGuess(fake_mep);
    for (size_t i = 0; i < size; ++i) {
      INFO("guess(" << i << ") = " << guess(i));
    }
    REQUIRE(guess.sum() == Approx(0.0));
  }
  SECTION("Uniform initial guess")
  {
    VIEFSolver solver(VPCMSolver::Uniform, VPCMSolver::SSD);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);
    Eigen::VectorXd guess = solver.initialGuess(fake_mep, total_charge);
    for (size_t i = 0; i < size; ++i) {
      INFO("guess(" << i << ") = " << guess(i));
    }
    REQUIRE(guess.sum()*nr_irrep == Approx(-total_charge));
  }
  SECTION("Diagonal initial guess")
  {
    VIEFSolver solver(VPCMSolver::Diagonal, VPCMSolver::SSD);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);
    Eigen::VectorXd guess = solver.initialGuess(fake_mep);
    for (size_t i = 0; i < size; ++i) {
      INFO("guess(" << i << ") = " << guess(i));
    }
    double reference = -170.9258887116;
    REQUIRE(guess.sum()*nr_irrep == Approx(reference));
  }
  SECTION("Low accuracy initial guess")
  {
    VIEFSolver solver(VPCMSolver::LowAccuracy, VPCMSolver::SSD);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);
    Eigen::VectorXd guess = solver.initialGuess(fake_mep);
    for (size_t i = 0; i < size; ++i) {
      INFO("guess(" << i << ") = " << guess(i));
    }
    double totalASC = - (2.0 * Ccharge + 4.0 * Hcharge) * (permittivity - 1) / permittivity;
    double guess_sum = guess.sum() * nr_irrep;
    CAPTURE(totalASC);
    CAPTURE(guess_sum);
    CAPTURE(totalASC - guess_sum);
    REQUIRE(guess_sum == Approx(totalASC).epsilon(1.0e-01));
  }
}
