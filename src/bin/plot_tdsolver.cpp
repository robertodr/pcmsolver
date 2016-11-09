/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
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

#include "bi_operators/Collocation.hpp"
#include "cavity/GePolCavity.hpp"
#include "OnsagerReactionField.hpp"
#include "td_solver/TDCPCMSolver.hpp"
#include "td_solver/TDIEFSolver.hpp"
#include "td_solver/TDSingleIEFSolver.hpp"
#include "green/Vacuum.hpp"
#include "utils/cnpy.hpp"
#include "utils/MathUtils.hpp"

extern "C" void host_writer(const char * message, size_t message_length);

void plot_tdief_lowdin_collocation(const GePolCavity &, const IGreensFunction &,
                                   const BoundaryIntegralOperator &, double, double,
                                   double, double, double, double);
void plot_tdcpcm_collocation(const GePolCavity &, const IGreensFunction &,
                             const BoundaryIntegralOperator &, double, double,
                             double, double, double);
void plot_tdsingleief_collocation(const GePolCavity &, const IGreensFunction &,
                                  const BoundaryIntegralOperator &, double, double,
                                  double, double, double);
void plot_tdonsagerief_collocation(const GePolCavity &, const IGreensFunction &,
                                   const BoundaryIntegralOperator &, double, double,
                                   double, double, double);

int main() {
  // global setup...
  initBohrToAngstrom(bohrToAngstrom);
  double radius = 1.181 * 1.10 / bohrToAngstrom();
  Sphere point(Eigen::Vector3d::Zero(), radius);
  double area = 0.4;
  double probeRadius = 0.0;
  double minRadius = 100.0;
  GePolCavity cavity(point, area, probeRadius, minRadius);

  Vacuum<> gfInside;
  integrator::Collocation biop;

  double e_0 = 35.69;
  double e_d = 1.807;
  double tau = 2000;                           // 48.38 fs
  double dt = 0.2;                             // 4.838 as
  double total_time = 100e-15 / secondsToAU(); // Total simulation time: 100 fs

  plot_tdief_lowdin_collocation(cavity, gfInside, biop, radius, e_0, e_d, tau, dt,
                                total_time);
  plot_tdcpcm_collocation(cavity, gfInside, biop, e_0, e_d, tau, dt, total_time);
  plot_tdsingleief_collocation(cavity, gfInside, biop, e_0, e_d, tau, dt,
                               total_time);
  plot_tdonsagerief_collocation(cavity, gfInside, biop, e_0, e_d, tau, dt,
                                total_time);

  return EXIT_SUCCESS;
}

void plot_tdief_lowdin_collocation(const GePolCavity & cavity,
                                   const IGreensFunction & gfInside,
                                   const BoundaryIntegralOperator & biop,
                                   double radius, double e_0, double e_d, double tau,
                                   double dt, double total_time) {
  int steps = int(total_time / dt) + 1; // Number of steps
  // The point-like dipole is at the origin, this is the direction
  Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
  bool cholesky = false;
  TDIEFSolver solver(e_0, e_d, tau, cholesky);
  solver.buildSystemMatrix(cavity, gfInside, biop);

  // Print diagonal matrices and their analytic form
  std::ofstream out;
  out.open("matrices_tdief_lowdin_collocation.dat");
  out.precision(16);
  for (size_t i = 0; i < cavity.size(); ++i) {
    out << solver.Lambda(i) << "  " << Lambda_lm(i) << "  " << solver.K_0(i) << "  "
        << K_lm(e_0, i) << "  " << solver.K_d(i) << "  " << K_lm(e_d, i) << "  "
        << solver.tau(i) * AUToFemtoseconds() << "  "
        << tau_lm(e_0, e_d, tau, i) * AUToFemtoseconds() << std::endl;
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
    out << "  " << t * secondsToAU() / 1.0e-15 << "  "
        << reactionField(radius, e_0, e_d, tau, t) << "  "
        << -fake_mep.dot(ASC_previous) << std::endl;
    ASC_previous = asc;
  }
  out.close();
}

void plot_tdcpcm_collocation(const GePolCavity & cavity,
                             const IGreensFunction & gfInside,
                             const BoundaryIntegralOperator & biop, double e_0,
                             double e_d, double tau, double dt, double total_time) {
  int steps = int(total_time / dt) + 1; // Number of steps
  // The point-like dipole is at the origin, this is the direction
  Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
  double corr = 0.0;
  TDCPCMSolver solver(e_0, e_d, tau, corr);
  solver.buildSystemMatrix(cavity, gfInside, biop);

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
    out << "  " << t * secondsToAU() / 1.0e-15 << "  " << -fake_mep.dot(ASC_previous)
        << std::endl;
    ASC_previous = asc;
  }
  out.close();
}

void plot_tdsingleief_collocation(const GePolCavity & cavity,
                                  const IGreensFunction & gfInside,
                                  const BoundaryIntegralOperator & biop, double e_0,
                                  double e_d, double tau, double dt,
                                  double total_time) {
  int steps = int(total_time / dt) + 1; // Number of steps
  // The point-like dipole is at the origin, this is the direction
  Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
  // Quadrupole relaxation time
  double tauIEF = tau * (3 * e_d + 2) / (3 * e_0 + 2);
  TDSingleIEFSolver solver(e_0, e_d, tau, tauIEF);
  solver.buildSystemMatrix(cavity, gfInside, biop);

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
    out << "  " << t * secondsToAU() / 1.0e-15 << "  " << -fake_mep.dot(ASC_previous)
        << std::endl;
    ASC_previous = asc;
  }
  out.close();
}

void plot_tdonsagerief_collocation(const GePolCavity & cavity,
                                   const IGreensFunction & gfInside,
                                   const BoundaryIntegralOperator & biop, double e_0,
                                   double e_d, double tau, double dt,
                                   double total_time) {
  int steps = int(total_time / dt) + 1; // Number of steps
  // The point-like dipole is at the origin, this is the direction
  Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
  double tauOnsager = tau * (2 * e_d + 1) / (2 * e_0 + 1);
  TDSingleIEFSolver solver(e_0, e_d, tau, tauOnsager);
  solver.buildSystemMatrix(cavity, gfInside, biop);

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
    out << "  " << t * secondsToAU() / 1.0e-15 << "  " << -fake_mep.dot(ASC_previous)
        << std::endl;
    ASC_previous = asc;
  }
  out.close();
}

void host_writer(const char * /* message */, size_t /* message_length */) {}
