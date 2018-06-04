/*
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

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

#include "Config.hpp"

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Collocation.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/Vacuum.hpp"
#include "td_solver/TDCPCMSolver.hpp"
#include "td_solver/TDIEFSolver.hpp"
#include "td_solver/TDSingleIEFSolver.hpp"
#include "utils/MathUtils.hpp"
#include "utils/cnpy.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using td_solver::TDCPCMSolver;
using td_solver::TDIEFSolver;
using td_solver::TDSingleIEFSolver;

void save_tdief_lowdin_collocation(const GePolCavity &,
                                   const IGreensFunction &,
                                   const Collocation &,
                                   double,
                                   double,
                                   double,
                                   double,
                                   double);
void save_tdcpcm_collocation(const GePolCavity &,
                             const IGreensFunction &,
                             const Collocation &,
                             double,
                             double,
                             double,
                             double,
                             double);
void save_tdsingleief_collocation(const GePolCavity &,
                                  const IGreensFunction &,
                                  const Collocation &,
                                  double,
                                  double,
                                  double,
                                  double,
                                  double);
void save_tdonsagerief_collocation(const GePolCavity &,
                                   const IGreensFunction &,
                                   const Collocation &,
                                   double,
                                   double,
                                   double,
                                   double,
                                   double);

int main() {
  initBohrToAngstrom(bohrToAngstrom);
  double radius = (1.181 * 1.10) / bohrToAngstrom();
  Molecule point = dummy<0>(radius);
  double area = 0.4;
  double probeRadius = 0.0;
  double minRadius = 100.0;
  GePolCavity cavity(point, area, probeRadius, minRadius);

  Vacuum<> gfInside;
  Collocation biop;

  double e_0 = 35.69;
  double e_d = 1.807;
  double tau = 2000;                           // 48.38 fs
  double dt = 0.2;                             // 4.838 as
  double total_time = 100e-15 / secondsToAU(); // Total simulation time: 100 fs

  save_tdief_lowdin_collocation(
      cavity, gfInside, biop, e_0, e_d, tau, dt, total_time);
  save_tdcpcm_collocation(cavity, gfInside, biop, e_0, e_d, tau, dt, total_time);
  save_tdsingleief_collocation(
      cavity, gfInside, biop, e_0, e_d, tau, dt, total_time);
  save_tdonsagerief_collocation(
      cavity, gfInside, biop, e_0, e_d, tau, dt, total_time);

  return EXIT_SUCCESS;
}

void save_tdief_lowdin_collocation(const GePolCavity & cavity,
                                   const IGreensFunction & gfInside,
                                   const Collocation & biop,
                                   double e_0,
                                   double e_d,
                                   double tau,
                                   double dt,
                                   double total_time) {
  int steps = int(total_time / dt) + 1; // Number of steps
  bool cholesky = false;
  TDIEFSolver solver(e_0, e_d, tau, cholesky);
  solver.buildSystemMatrix(cavity, gfInside, biop);

  // Save diagonal matrix entries
  int size = cavity.size();
  unsigned int dim = static_cast<unsigned int>(size);
  cnpy::custom::npy_save("tdief_lowdin_collocation_Lambda.npy", solver.Lambda());
  cnpy::custom::npy_save("tdief_lowdin_collocation_K_0.npy", solver.K_0());
  cnpy::custom::npy_save("tdief_lowdin_collocation_K_d.npy", solver.K_d());
  cnpy::custom::npy_save("tdief_lowdin_collocation_tau.npy", solver.tau());

  // The point-like dipole is at the origin, this is the direction
  Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
  Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < size; ++i) {
    Eigen::Vector3d center = cavity.elementCenter(i);
    double distance = center.norm();
    fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
  }
  // Time-dependent IEF reaction field
  std::vector<double> TDIEF;
  TDIEF.reserve(steps);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
  for (int i = 0; i < steps; ++i) {
    TDIEF.push_back(-fake_mep.dot(ASC_previous));
    asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
    ASC_previous = asc;
  }

  dim = static_cast<unsigned int>(steps);
  const unsigned int new_shape[] = {dim};
  cnpy::npy_save(
      "tdief_lowdin_collocation.npy", TDIEF.data(), new_shape, 1, "w", true);
}

void save_tdcpcm_collocation(const GePolCavity & cavity,
                             const IGreensFunction & gfInside,
                             const Collocation & biop,
                             double e_0,
                             double e_d,
                             double tau,
                             double dt,
                             double total_time) {
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
  std::vector<double> TDCPCM;
  TDCPCM.reserve(steps);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
  for (int i = 0; i < steps; ++i) {
    TDCPCM.push_back(-fake_mep.dot(ASC_previous));
    asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
    ASC_previous = asc;
  }

  unsigned int dim = static_cast<unsigned int>(steps);
  const unsigned int shape[] = {dim};
  cnpy::npy_save("tdcpcm_collocation.npy", TDCPCM.data(), shape, 1, "w", true);
}

void save_tdsingleief_collocation(const GePolCavity & cavity,
                                  const IGreensFunction & gfInside,
                                  const Collocation & biop,
                                  double e_0,
                                  double e_d,
                                  double tau,
                                  double dt,
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
  std::vector<double> TDIEF;
  TDIEF.reserve(steps);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
  for (int i = 0; i < steps; ++i) {
    TDIEF.push_back(-fake_mep.dot(ASC_previous));
    asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
    ASC_previous = asc;
  }

  unsigned int dim = static_cast<unsigned int>(steps);
  const unsigned int shape[] = {dim};
  cnpy::npy_save("tdsingleief_collocation.npy", TDIEF.data(), shape, 1, "w", true);
}

void save_tdonsagerief_collocation(const GePolCavity & cavity,
                                   const IGreensFunction & gfInside,
                                   const Collocation & biop,
                                   double e_0,
                                   double e_d,
                                   double tau,
                                   double dt,
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
  std::vector<double> TDIEF;
  TDIEF.reserve(steps);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
  for (int i = 0; i < steps; ++i) {
    TDIEF.push_back(-fake_mep.dot(ASC_previous));
    asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
    ASC_previous = asc;
  }

  unsigned int dim = static_cast<unsigned int>(steps);
  const unsigned int shape[] = {dim};
  cnpy::npy_save("tdonsagerief_collocation.npy", TDIEF.data(), shape, 1, "w", true);
}
