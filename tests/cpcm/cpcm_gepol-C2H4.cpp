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

#include <catch.hpp>

#include <cmath>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "CPCMSolver.hpp"
#include "TestingMolecules.hpp"

SCENARIO("Test solver for the C-PCM and the C2H4 molecule", "[solver][cpcm][cpcm_symmetry][cpcm_gepol-C2H4]")
{
    GIVEN("An isotropic environment modelled as a perfect conductor and a point charge")
    {
        double permittivity = 78.39;
        Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
        UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
            UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
        bool symm = true;
        double correction = 0.0;

        double Ccharge = 6.0;
        double Hcharge = 1.0;
        // The total ASC for a conductor is -Q
        // for CPCM it will be -Q*(epsilon-1)/epsilon
        double totalASC = - (2.0 * Ccharge + 4.0 * Hcharge) * (permittivity - 1) / permittivity;

        double area = 1.5;
        double probeRadius = 1.385;
        double minRadius = 100.0;

        /*! \class CPCMSolver
         *  \test \b C2H4GePolC1 tests CPCMSolver using C2H4 with a GePol cavity in C1 symmetry
         */
        WHEN("the point group is C1")
        {
          Molecule molec = C2H4_C1();
          GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius, "cpcm_c1_noadd");

          CPCMSolver solver(symm, correction);
          solver.buildSystemMatrix(cavity, gfInside, gfOutside);

          size_t size = cavity.size();
          Eigen::VectorXd mep_C1 = computeMEP(molec, cavity.elements());
          for (size_t i = 0; i < size; ++i) {
            INFO("mep_C1(" << i << ") = " << mep_C1(i));
          }

          THEN("the total apparent surface charge is")
          {
            Eigen::VectorXd asc_C1 = solver.computeCharge(mep_C1);

            for (size_t i = 0; i < size; ++i) {
              INFO("asc_C1(" << i << ") = " << asc_C1(i));
            }

            int nr_irrep = cavity.pointGroup().nrIrrep();
            double totalFakeASC = asc_C1.sum() * nr_irrep;
            CAPTURE(totalASC);
            CAPTURE(totalFakeASC);
            CAPTURE(totalASC - totalFakeASC);
            REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
          }
        }
        /*! \class CPCMSolver
         *  \test \b C2H4GePolD2h tests CPCMSolver using C2H4 with a GePol cavity in D2h symmetry
         */
        WHEN("the point group is D2h")
        {
          Molecule molec = C2H4_D2h();
          GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius, "cpcm_d2h_noadd");

          CPCMSolver solver(symm, correction);
          solver.buildSystemMatrix(cavity, gfInside, gfOutside);

          size_t size = cavity.size();
          Eigen::VectorXd mep_D2h = computeMEP(molec, cavity.elements());
          for (size_t i = 0; i < size; ++i) {
            INFO("mep_D2h(" << i << ") = " << mep_D2h(i));
          }

          THEN("the total apparent surface charge is")
          {
            Eigen::VectorXd asc_D2h = solver.computeCharge(mep_D2h);

            for (size_t i = 0; i < size; ++i) {
              INFO("asc_D2h(" << i << ") = " << asc_D2h(i));
            }

            int nr_irrep = cavity.pointGroup().nrIrrep();
            double totalFakeASC = asc_D2h.sum() * nr_irrep;
            CAPTURE(totalASC);
            CAPTURE(totalFakeASC);
            CAPTURE(totalASC - totalFakeASC);
            REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
          }
        }
    }
}
