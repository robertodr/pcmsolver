/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "bi_operators/BIOperatorData.hpp"
#include "cavity/CavityData.hpp"
#include "green/GreenData.hpp"
#include "interface/Input.hpp"
#include "solver/SolverData.hpp"
#include "td_solver/TDSolverData.hpp"
#include "utils/Sphere.hpp"

using pcm::Input;
using pcm::utils::Sphere;

/*! \class Input
 *  \test \b InputTDSolverTest_TDSolver tests input reading on an input file parsed
 * by pcmsolver.py
 */
TEST_CASE("Input reading using GetKw for an input file for a TD solver with solvent "
          "chosen by name",
          "[input][tdsolver][input_tdsolver]") {
  std::string filename = "@tdsolver.inp";
  Input parsedInput(filename);
  std::string units = "ANGSTROM";
  int CODATAyear = 2002;
  std::string type = "GEPOL";
  std::string radiiSet = "BONDI";
  std::string mode = "IMPLICIT";
  std::string solverType = "IEFPCM";
  std::string TDsolverType = "TDSINGLEIEF";
  std::string greenInsideType = "VACUUM_DERIVATIVE";
  std::string greenOutsideType = "UNIFORMDIELECTRIC_DERIVATIVE";
  bool initWithDynamic = false;
  int derivativeInsideType = 1;
  int derivativeOutsideType = 1;
  double area = 10.0 * angstrom2ToBohr2();
  std::string solvent = "Acetonitrile";     // Name in the Solvent object
  double tauIEF = 150 / 2.418884326505e-02; // in au
  double tau = 48.38 / 2.418884326505e-02;  // in au
  bool isTD = true;
  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(type == parsedInput.cavityParams().cavityType);
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(solverType == parsedInput.solverParams().solverType);
  REQUIRE(TDsolverType == parsedInput.TDSolverParams().TDsolverType);
  REQUIRE(greenInsideType == parsedInput.insideGreenParams().greensFunctionType);
  REQUIRE(greenOutsideType ==
          parsedInput.outsideStaticGreenParams().greensFunctionType);
  REQUIRE(initWithDynamic == parsedInput.TDSolverParams().initWithDynamic);
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(solvent == parsedInput.solvent().name);
  REQUIRE(tauIEF == Approx(parsedInput.TDSolverParams().tauIEF));
  REQUIRE(tau == Approx(parsedInput.TDSolverParams().tau));
  REQUIRE(isTD == parsedInput.isTD());
}

/*! \class Input
 *  \test \b InputExplicitTDSolverTest_TDSolver tests input reading on an input file
 * parsed by pcmsolver.py
 */
TEST_CASE("Input reading using GetKw for an input file for a TD solver with solvent "
          "explicitly given",
          "[input][tdsolver][explicit][input_tdsolver]") {
  std::string filename = "@explicit_td.inp";
  Input parsedInput(filename);
  std::string units = "ANGSTROM";
  int CODATAyear = 2002;
  std::string type = "GEPOL";
  std::string radiiSet = "BONDI";
  std::string mode = "IMPLICIT";
  std::string solverType = "IEFPCM";
  std::string TDsolverType = "TDSINGLEIEF";
  std::string greenInsideType = "VACUUM_DERIVATIVE";
  std::string greenOutsideType = "UNIFORMDIELECTRIC_DERIVATIVE";
  bool initWithDynamic = true;
  int derivativeInsideType = 1;
  int derivativeOutsideType = 1;
  double area = 10.0 * angstrom2ToBohr2();
  std::string solvent = "Acetonitrile";     // Name in the Solvent object
  double tauIEF = 150 / 2.418884326505e-02; // in au
  double tau = 48.38 / 2.418884326505e-02;  // in au
  double epsilonInside = 1.0;
  double epsilonStaticOutside = 78.39;
  double epsilonDynamicOutside = 2.230;
  bool isTD = true;
  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(type == parsedInput.cavityParams().cavityType);
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(solverType == parsedInput.solverParams().solverType);
  REQUIRE(TDsolverType == parsedInput.TDSolverParams().TDsolverType);
  REQUIRE(greenInsideType == parsedInput.insideGreenParams().greensFunctionType);
  REQUIRE(greenOutsideType ==
          parsedInput.outsideStaticGreenParams().greensFunctionType);
  REQUIRE(initWithDynamic == parsedInput.TDSolverParams().initWithDynamic);
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(tauIEF == Approx(parsedInput.TDSolverParams().tauIEF));
  REQUIRE(tau == Approx(parsedInput.TDSolverParams().tau));
  REQUIRE(epsilonInside == Approx(parsedInput.insideGreenParams().epsilon));
  REQUIRE(epsilonStaticOutside ==
          Approx(parsedInput.outsideStaticGreenParams().epsilon));
  REQUIRE(epsilonDynamicOutside ==
          Approx(parsedInput.outsideDynamicGreenParams().epsilon));
  REQUIRE(isTD == parsedInput.isTD());
}
