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
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "VIEFSolver.hpp"
#include "TestingMolecules.hpp"

SCENARIO("Test variational solver for the anisotropic IEFPCM for a point charge in different Abelian point groups", "[variational_solver][viefpcm][viefpcm_anisotropic-symmetry][anisotropic]")
{
    GIVEN("An isotropic environment modelled and a point charge forcing the use of an anisotropic solver")
    {
        double permittivity = 78.39;
        Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
        UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
            UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
        bool symm = true;

        double charge = 8.0;
        double totalASC = - charge * (permittivity - 1) / permittivity;

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolC1 tests VIEFSolver using a point charge with a GePol cavity in C1 symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is C1")
        {
            Molecule point = dummy<0>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "C1");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolC2 tests VIEFSolver using a point charge with a GePol cavity in C2 symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is C2")
        {
            Molecule point = dummy<1>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "C2");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolCs tests VIEFSolver using a point charge with a GePol cavity in Cs symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is Cs")
        {
            Molecule point = dummy<2>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "Cs");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolCi tests VIEFSolver using a point charge with a GePol cavity in Ci symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is Ci")
        {
            Molecule point = dummy<3>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "Ci");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolD2 tests VIEFSolver using a point charge with a GePol cavity in D2 symmetry
         */
        WHEN("the point group is D2")
        {
            Molecule point = dummy<4>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "D2");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolC2v tests VIEFSolver using a point charge with a GePol cavity in C2v symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is C2v")
        {
            Molecule point = dummy<5>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "C2v");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolC2h tests VIEFSolver using a point charge with a GePol cavity in C2h symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is C2h")
        {
            Molecule point = dummy<6>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "C2h");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }

        /*! \class VIEFSolver
         *  \test \b pointChargeGePolD2h tests VIEFSolver using a point charge with a GePol cavity in D2h symmetry
         *  We are forcing the usage of the buildAnisotropicMatrix method.
         *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
         *  The point charge is at the origin.
         */
        WHEN("the point group is D2h")
        {
            Molecule point = dummy<7>(2.929075493);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius, "D2h");

            VIEFSolver aniso_solver;
            aniso_solver.buildAnisotropicMatrix(cavity, gfInside, gfOutside);

            VIEFSolver iso_solver;
            iso_solver.buildIsotropicMatrix(cavity, gfInside, gfOutside);

            size_t size = cavity.size();
            Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

            THEN("the total apparent surface charge is")
            {
                Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
                Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

                int nr_irrep = cavity.pointGroup().nrIrrep();
                double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
                double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

                CAPTURE(totalASC);
                CAPTURE(totalAnisoASC);
                CAPTURE(totalASC - totalAnisoASC);
                REQUIRE(totalIsoASC == Approx(totalAnisoASC));
                REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
            }
        }
    }
}
