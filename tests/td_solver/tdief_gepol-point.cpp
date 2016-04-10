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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>

#include "Config.hpp"

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/CollocationIntegrator.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/Vacuum.hpp"
#include "green/UniformDielectric.hpp"
#include "solver/IEFSolver.hpp"
#include "td_solver/TDIEFSolver.hpp"
#include "utils/cnpy.hpp"
#include "utils/MathUtils.hpp"

SCENARIO("Test time-dependent solver for the IEFPCM for a point dipole and a GePol cavity", "[td_solver][tdiefpcm][tdief_gepol-point]")
{
    GIVEN("An isotropic environment and a point dipole")
    {
        double e_0 = 35.69;
        double e_d = 1.807;
        double tau = 2000; // 48.38 fs
        double dt = 0.2; // 4.838 as
        double secondsToAU = 2.418884326509e-17;
        double total_time = 100e-15 / secondsToAU; // Total simulation time: 100 fs
        int steps = int(total_time/dt) + 1; // Number of steps
        bool cholesky = false;

        /*! \class TDIEFSolver
         *  \test \b pointDipoleGePolLowdin tests TDIEFSolver using a point dipole with a GePol cavity
         *  \note uses Lowdin diagonalization procedure
         */
        WHEN("the point dipole is located at the origin")
        {
            double radius = 1.181 * 1.10 / bohrToAngstrom();
            Molecule point = dummy<0>(radius);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius);
            size_t size = cavity.size();

            TDIEFSolver<> solver(e_0, e_d, tau, cholesky);
            solver.buildSystemMatrix(cavity);

            THEN("the elements of the Lambda matrix are")
            {
                Eigen::VectorXd Lambda_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_Lambda.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.Lambda(i) == Approx(Lambda_ref(i)));
                }
            }
            AND_THEN("the elements of the K_0 matrix are")
            {
                Eigen::VectorXd K_0_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_K_0.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.K_0(i) == Approx(K_0_ref(i)));
                }
            }
            AND_THEN("the elements of the K_d matrix are")
            {
                Eigen::VectorXd K_d_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_K_d.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.K_d(i) == Approx(K_d_ref(i)));
                }
            }
            AND_THEN("the elements of the tau matrix are")
            {
                Eigen::VectorXd tau_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_tau.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.tau(i) == Approx(tau_ref(i)));
                }
            }

            Vacuum<> gfInside;
            UniformDielectric<> gfOutside(e_d);
            bool symm = true;
            IEFSolver static_solver(symm);
            static_solver.buildSystemMatrix(cavity, gfInside, gfOutside);

            // The point-like dipole is at the origin, this is the direction
            Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
            Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
            for (size_t i = 0; i < size; ++i) {
                Eigen::Vector3d center = cavity.elementCenter(i);
                double distance = center.norm();
                fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
            }
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }
            Eigen::VectorXd fake_asc = static_solver.computeCharge(fake_mep);

            Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
            for (size_t i = 0; i < size; ++i) {
                INFO("ASC_previous(" << i << ") = " << ASC_previous(i));
            }
            // Internal consistency check
            AND_THEN("the ASC at t = 0 is the dynamic ASC")
            {
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(fake_asc(i) == Approx(ASC_previous(i)));
                }
            }

            // Time-dependent IEF reaction field
            std::vector<double> TDIEF;
            TDIEF.reserve(steps);
            Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
            for (int i = 0; i < steps; ++i) {
                TDIEF.push_back(-fake_mep.dot(ASC_previous));
                asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
                ASC_previous = asc;
            }
            AND_THEN("the ASC trajectory is")
            {
                cnpy::NpyArray raw_TDIEF_ref = cnpy::npy_load("tdief_lowdin_collocation.npy");
                double * TDIEF_ref = reinterpret_cast<double *>(raw_TDIEF_ref.data);
                for (int i = 0; i < steps; ++i) {
                    REQUIRE(TDIEF[i] == Approx(TDIEF_ref[i]));
                }
            }
        }

        /*! \class TDIEFSolver
         *  \test \b pointDipoleGePolLowdin tests TDIEFSolver using a point dipole with a GePol cavity
         *  The point dipole is away from the origin.
         *  \note uses Lowdin diagonalization procedure
         */
        AND_WHEN("the point dipole is located away from the origin")
        {
            Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
            double radius = 1.181 * 1.10 / bohrToAngstrom();
            Molecule point = dummy<0>(radius, origin);
            double area = 0.4;
            double probeRadius = 0.0;
            double minRadius = 100.0;
            GePolCavity cavity(point, area, probeRadius, minRadius);
            size_t size = cavity.size();

            TDIEFSolver<> solver(e_0, e_d, tau, cholesky);
            solver.buildSystemMatrix(cavity);

            THEN("the elements of the Lambda matrix are")
            {
                Eigen::VectorXd Lambda_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_Lambda.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.Lambda(i) == Approx(Lambda_ref(i)));
                }
            }
            AND_THEN("the elements of the K_0 matrix are")
            {
                Eigen::VectorXd K_0_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_K_0.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.K_0(i) == Approx(K_0_ref(i)));
                }
            }
            AND_THEN("the elements of the K_d matrix are")
            {
                Eigen::VectorXd K_d_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_K_d.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.K_d(i) == Approx(K_d_ref(i)));
                }
            }
            AND_THEN("the elements of the tau matrix are")
            {
                Eigen::VectorXd tau_ref = cnpy::custom::npy_load<double>("tdief_lowdin_collocation_tau.npy");
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(solver.tau(i) == Approx(tau_ref(i)));
                }
            }

            Vacuum<> gfInside;
            UniformDielectric<> gfOutside(e_d);
            bool symm = true;
            IEFSolver static_solver(symm);
            static_solver.buildSystemMatrix(cavity, gfInside, gfOutside);

            // The point-like dipole is at the origin, this is the direction
            Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
            Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
            for (size_t i = 0; i < size; ++i) {
                Eigen::Vector3d center = cavity.elementCenter(i);
                double distance = (center - origin).norm();
                fake_mep(i) = (center - origin).dot(dipole) / std::pow(distance, 3);
            }
            for (size_t i = 0; i < size; ++i) {
                INFO("fake_mep(" << i << ") = " << fake_mep(i));
            }
            Eigen::VectorXd fake_asc = static_solver.computeCharge(fake_mep);

            Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
            for (size_t i = 0; i < size; ++i) {
                INFO("ASC_previous(" << i << ") = " << ASC_previous(i));
            }
            // Internal consistency check
            AND_THEN("the ASC at t = 0 is the dynamic ASC")
            {
                for (size_t i = 0; i < size; ++i) {
                    REQUIRE(fake_asc(i) == Approx(ASC_previous(i)));
                }
            }

            // Time-dependent IEF reaction field
            std::vector<double> TDIEF;
            TDIEF.reserve(steps);
            Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
            for (int i = 0; i < steps; ++i) {
                TDIEF.push_back(-fake_mep.dot(ASC_previous));
                asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
                ASC_previous = asc;
            }
            AND_THEN("the ASC trajectory is")
            {
                cnpy::NpyArray raw_TDIEF_ref = cnpy::npy_load("tdief_lowdin_collocation.npy");
                double * TDIEF_ref = reinterpret_cast<double *>(raw_TDIEF_ref.data);
                for (int i = 0; i < steps; ++i) {
                    REQUIRE(TDIEF[i] == Approx(TDIEF_ref[i]));
                }
            }
        }
    }
}
