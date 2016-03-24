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

#include <vector>
#include <cmath>


#include <Eigen/Core>

#include "cavity/GePolCavity.hpp"
#include "bi_operators/CollocationIntegrator.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/Vacuum.hpp"
#include "green/IonicLiquid.hpp"
#include "green/UniformDielectric.hpp"
#include "utils/ChargeDistribution.hpp"
#include "TestingMolecules.hpp"

SCENARIO("Calculation of the Newton potential", "[utils][newton_potential][utils_newton-potential]")
{
  GIVEN("A classical charge distribution of monopoles")
  {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    Molecule molec = C2H4();
    GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius, "newton");
    ChargeDistribution dist(molec.charges(), molec.geometry());

    WHEN("the Newton potential is calculated in vacuum")
    {
      Vacuum<AD_directional, CollocationIntegrator> gf = Vacuum<AD_directional, CollocationIntegrator>();
      Eigen::VectorXd newton = computeNewtonPotential(
          pcm::bind(&Vacuum<AD_directional, CollocationIntegrator>::kernelS, gf, pcm::_1, pcm::_2),
          cavity.elementCenter(), dist);
      THEN("comparison with the computeMEP method results in")
      {
        Eigen::VectorXd mep = computeMEP(molec, cavity.elements());
        REQUIRE(newton.size() == mep.size());
        for (int i = 0; i < mep.size(); ++i) {
          REQUIRE(newton(i) == Approx(mep(i)));
        }
      }
    }

    AND_WHEN("the Newton potential is calculated in a uniform dielectric")
    {
      double permittivity = 78.39;
      UniformDielectric<AD_directional, CollocationIntegrator> gf =
        UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
      Eigen::VectorXd newton = computeNewtonPotential(
          pcm::bind(&UniformDielectric<AD_directional, CollocationIntegrator>::kernelS, gf, pcm::_1, pcm::_2),
          cavity.elementCenter(), dist);
      THEN("comparison with the computeMEP method results in")
      {
        Eigen::VectorXd mep = computeMEP(molec, cavity.elements());
        mep /= permittivity;
        REQUIRE(newton.size() == mep.size());
        for (int i = 0; i < mep.size(); ++i) {
          REQUIRE(newton(i) == Approx(mep(i)));
        }
      }
    }

    AND_WHEN("the Newton potential is calculated in an ionic liquid with kappa = 0.0")
    {
      double permittivity = 78.39;
      double kappa = 0.0;
      IonicLiquid<AD_directional, CollocationIntegrator> gf =
        IonicLiquid<AD_directional, CollocationIntegrator>(permittivity, kappa);
      Eigen::VectorXd newton = computeNewtonPotential(
          pcm::bind(&IonicLiquid<AD_directional, CollocationIntegrator>::kernelS, gf, pcm::_1, pcm::_2),
          cavity.elementCenter(), dist);
      THEN("comparison with the computeMEP method results in")
      {
        Eigen::VectorXd mep = computeMEP(molec, cavity.elements());
        mep /= permittivity;
        REQUIRE(newton.size() == mep.size());
        for (int i = 0; i < mep.size(); ++i) {
          REQUIRE(newton(i) == Approx(mep(i)));
        }
      }
    }
  }
}
