/**
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

#include "ddPCM.hpp"

#include <cmath>

#include "utils/Molecule.hpp"

namespace pcm {
namespace solver {
ddPCM::ddPCM(const Molecule & m) : nSpheres_(m.spheres().size()), molecule_(m) {
  int ncav = 0;
  nSpheres_ = m.spheres().size();
  double * xs = new double[nSpheres_];
  double * ys = new double[nSpheres_];
  double * zs = new double[nSpheres_];
  double * rs = new double[nSpheres_];
  for (int i = 0; i < nSpheres_; ++i) {
    xs[i] = m.spheres(i).center(0);
    ys[i] = m.spheres(i).center(1);
    zs[i] = m.spheres(i).center(2);
    rs[i] = m.spheres(i).radius;
  }
  ddinit(&nSpheres_, xs, ys, zs, rs, &ncav);
  cavity_ = Eigen::Matrix3Xd::Zero(3, ncav);
  copy_cavity(cavity_.data());
  delete[] xs;
  delete[] ys;
  delete[] zs;
  delete[] rs;
}

ddPCM::~ddPCM() { memfree(); }

Eigen::MatrixXd ddPCM::computeCharges(const Eigen::VectorXd & phi) const {
  int nbasis = 7 * 7;
  double ene = 0.0;
  Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(nbasis, nSpheres_);
  Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(nbasis, nSpheres_);
  for (int i = 0; i < nSpheres_; ++i) {
    psi(0, i) = std::sqrt(4.0 * M_PI) * molecule_.charges(i);
  }
  itsolv_direct(phi.data(), psi.data(), sigma.data(), &ene);
  return sigma;
}
} // namespace solver
} // namespace pcm
