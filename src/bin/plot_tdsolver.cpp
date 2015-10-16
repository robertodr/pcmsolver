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

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Core>

#include "cnpy.hpp"
#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "OnsagerReactionField.hpp"
#include "PhysicalConstants.hpp"
#include "TDCPCMSolver.hpp"
#include "TDIEFSolver.hpp"
#include "TDSingleIEFSolver.hpp"
#include "TDOnsagerIEFSolver.hpp"
#include "TestingMolecules.hpp"
#include "Vacuum.hpp"

const double secondsToAU = 2.418884326509e-17;
const double AUTofemtoseconds = secondsToAU / 1.0e-15;

void plot_tdief_lowdin_collocation(const GePolCavity &, double, double, double, double, double, double);
void plot_tdcpcm_collocation(const GePolCavity &, double, double, double, double, double);
void plot_tdsingleief_collocation(const GePolCavity &, double, double, double, double, double);
void plot_tdonsagerief_collocation(const GePolCavity &, double, double, double, double, double);

int main()
{
    double radius = 1.181 * 1.10 / convertBohrToAngstrom;
    Molecule point = dummy<0>(radius);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);

    double e_0 = 35.69;
    double e_d = 1.807;
    double tau = 2000; // 48.38 fs
    double dt = 0.2; // 4.838 as
    double total_time = 100e-15 / secondsToAU; // Total simulation time: 100 fs

    plot_tdief_lowdin_collocation(cavity, radius, e_0, e_d, tau, dt, total_time);
    plot_tdcpcm_collocation(cavity, e_0, e_d, tau, dt, total_time);
    plot_tdsingleief_collocation(cavity, e_0, e_d, tau, dt, total_time);
    plot_tdonsagerief_collocation(cavity, e_0, e_d, tau, dt, total_time);

    return EXIT_SUCCESS;
}

void plot_tdief_lowdin_collocation(const GePolCavity & cavity, double radius, double e_0, double e_d, double tau, double dt, double total_time)
{
    int steps = int(total_time/dt) + 1; // Number of steps
    // The point-like dipole is at the origin, this is the direction
    Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
    bool cholesky = false;
    TDIEFSolver<AD_directional, CollocationIntegrator> solver(e_0, e_d, tau, cholesky);
    solver.buildSystemMatrix(cavity);

    // Print diagonal matrices and their analytic form
    std::ofstream out;
    out.open("matrices_tdief_lowdin_collocation.dat");
    out.precision(16);
    for (size_t i = 0; i < cavity.size(); ++i) {
        out << solver.Lambda(i) << "  "
            << Lambda_lm(i)  << "  "
            << solver.K_0(i) << "  "
            << K_lm(e_0, i)  << "  "
            << solver.K_d(i) << "  "
            << K_lm(e_d, i)  << "  "
            << solver.tau(i) * AUTofemtoseconds << "  "
            << tau_lm(e_0, e_d, tau, i) * AUTofemtoseconds << std::endl;
    }
    out.close();

    size_t size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (size_t i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
    }
    // Propagate
    double t_0 = 0.0, t = 0.0;
    Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
    out.open("tdief_lowdin_collocation.dat");
    out.precision(16);
    for (int i = 0; i < steps; ++i) {
        t = t_0 + i * dt;
        asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
        out << "  " << t * secondsToAU / 1.0e-15
            << "  " << reactionField(radius, e_0, e_d, tau, t)
            << "  " << -fake_mep.dot(ASC_previous)
            << std::endl;
        ASC_previous = asc;
    }
    out.close();
}

void plot_tdcpcm_collocation(const GePolCavity & cavity, double e_0, double e_d, double tau, double dt, double total_time)
{
    int steps = int(total_time/dt) + 1; // Number of steps
    // The point-like dipole is at the origin, this is the direction
    Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
    double corr = 0.0;
    TDCPCMSolver<AD_directional, CollocationIntegrator> solver(e_0, e_d, tau, corr);
    solver.buildSystemMatrix(cavity);

    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
    }
    // Propagate
    double t_0 = 0.0, t = 0.0;
    Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
    std::ofstream out;
    out.open("tdcpcm_collocation.dat");
    out.precision(16);
    for (int i = 0; i < steps; ++i) {
        t = t_0 + i * dt;
        asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
        out << "  " << t * secondsToAU / 1.0e-15
            << "  " << -fake_mep.dot(ASC_previous)
            << std::endl;
        ASC_previous = asc;
    }
    out.close();
}

void plot_tdsingleief_collocation(const GePolCavity & cavity, double e_0, double e_d, double tau, double dt, double total_time)
{
    int steps = int(total_time/dt) + 1; // Number of steps
    // The point-like dipole is at the origin, this is the direction
    Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
    // Quadrupole relaxation time
    double tauIEF = tau * (3 * e_d + 2) / (3 * e_0 + 2);
    TDSingleIEFSolver<AD_directional, CollocationIntegrator> solver(e_0, e_d, tau, tauIEF);
    solver.buildSystemMatrix(cavity);

    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
    }
    // Propagate
    double t_0 = 0.0, t = 0.0;
    Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
    std::ofstream out;
    out.open("tdsingleief_collocation.dat");
    out.precision(16);
    for (int i = 0; i < steps; ++i) {
        t = t_0 + i * dt;
        asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
        out << "  " << t * secondsToAU / 1.0e-15
            << "  " << -fake_mep.dot(ASC_previous)
            << std::endl;
        ASC_previous = asc;
    }
    out.close();
}

void plot_tdonsagerief_collocation(const GePolCavity & cavity, double e_0, double e_d, double tau, double dt, double total_time)
{
    int steps = int(total_time/dt) + 1; // Number of steps
    // The point-like dipole is at the origin, this is the direction
    Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
    TDOnsagerIEFSolver<AD_directional, CollocationIntegrator> solver(e_0, e_d, tau);
    solver.buildSystemMatrix(cavity);

    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
    }
    // Propagate
    double t_0 = 0.0, t = 0.0;
    Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
    std::ofstream out;
    out.open("tdonsagerief_collocation.dat");
    out.precision(16);
    for (int i = 0; i < steps; ++i) {
        t = t_0 + i * dt;
        asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
        out << "  " << t * secondsToAU / 1.0e-15
            << "  " << -fake_mep.dot(ASC_previous)
            << std::endl;
        ASC_previous = asc;
    }
    out.close();
}
