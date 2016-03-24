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


#include <Eigen/Core>

#include "bi_operators/CollocationIntegrator.hpp"
#include "green/DerivativeTypes.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/Vacuum.hpp"
#include "green/UniformDielectric.hpp"
#include "solver/IEFSolver.hpp"
#include "TestingMolecules.hpp"

SCENARIO("Test a point charge and a GePol cavity in flipped environment (uniform dielectric inside, vacuum outside)", "[solver][iefpcm][iefpcm_gepol-point_flipped][flipped]")
{
    GIVEN("An isotropic environment inside the cavity and a point charge")
    {
        double permittivity = 78.39;
        UniformDielectric<AD_directional, CollocationIntegrator> gf_i =
            UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
        Vacuum<AD_directional, CollocationIntegrator> gf_o = Vacuum<AD_directional, CollocationIntegrator>();
        bool symm = true;

        double charge = 8.0;
        double totalASC = charge * (permittivity - 1) / permittivity;

        /*! \class IEFSolver
         *  \test \b pointChargeGePolFlipped tests IEFSolver using a point charge with a GePol cavity and a flipped environment
         *  We force also usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point charge is located at the origin")
        {
            Molecule point = dummy<0>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);

            IEFSolver aniso_solver(symm);
            aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o);

            IEFSolver iso_solver(symm);
            iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }

            THEN("the apparent surface charge is")
            {
              // The RHS has the MEP scaled by the permittivity
              Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep / permittivity);
              // The RHS has the opposite sign for the flipped case
              Eigen::VectorXd iso_fake_asc = - iso_solver.computeCharge(fake_mep);
              double totalAnisoASC = aniso_fake_asc.sum();
              double totalIsoASC = iso_fake_asc.sum();

              for (size_t i = 0; i < size; ++i) {
                INFO("aniso_fake_asc(" << i << ") = " << aniso_fake_asc(i));
              }
              for (size_t i = 0; i < size; ++i) {
                INFO("iso_fake_asc(" << i << ") = " << iso_fake_asc(i));
              }
              print_eigen_matrix(aniso_fake_asc, "aniso.log");
              print_eigen_matrix(iso_fake_asc, "iso.log");

              CAPTURE(totalASC);
              CAPTURE(totalAnisoASC);
              CAPTURE(totalASC - totalAnisoASC);
              REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
              REQUIRE(totalIsoASC == Approx(totalAnisoASC));
            }
        }

        /*! \class IEFSolver
         *  \test \b pointChargeShiftedGePolFlipped tests IEFSolver using a point charge with a GePol cavity and a flipped environment
         *  We force also usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is away from the origin.
         */
        AND_WHEN("the point charge is located away from the origin")
        {
            Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
            Molecule point = dummy<0>(2.929075493, origin);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);

            IEFSolver aniso_solver(symm);
            aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o);

            IEFSolver iso_solver(symm);
            iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o);

            double charge = 8.0;
            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge, origin);
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }

            THEN("the apparent surface charge is")
            {
              // The RHS has the MEP scaled by the permittivity
              Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep / permittivity);
              // The RHS has the opposite sign for the flipped case
              Eigen::VectorXd iso_fake_asc = - iso_solver.computeCharge(fake_mep);
              double totalAnisoASC = aniso_fake_asc.sum();
              double totalIsoASC = iso_fake_asc.sum();

              for (size_t i = 0; i < size; ++i) {
                INFO("aniso_fake_asc(" << i << ") = " << aniso_fake_asc(i));
              }
              for (size_t i = 0; i < size; ++i) {
                INFO("iso_fake_asc(" << i << ") = " << iso_fake_asc(i));
              }

              CAPTURE(totalASC);
              CAPTURE(totalAnisoASC);
              CAPTURE(totalASC - totalAnisoASC);
              REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
              REQUIRE(totalIsoASC == Approx(totalAnisoASC));
            }
        }
    }
}
