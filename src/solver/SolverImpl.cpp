/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators.
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

#include "SolverImpl.hpp"

#include <cmath>

#include "Config.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>

#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"

namespace solver {
Eigen::MatrixXd computeS(const Cavity & cav, const IGreensFunction & gf) {
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute S
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    // Fill diagonal
    Element source = cav.elements(i);
    S(i, i) = gf.singleLayer(source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      // Fill off-diagonal
      Element probe = cav.elements(j);
      if (i != j)
        S(i, j) = gf.kernelS(source.center(), probe.center());
    }
  }

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(S, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }

  return S;
}

Eigen::MatrixXd computeD(const Cavity & cav, const IGreensFunction & gf) {
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute D
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    // Fill diagonal
    Element source = cav.elements(i);
    D(i, i) = gf.doubleLayer(source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      // Fill off-diagonal
      Element probe = cav.elements(j);
      if (i != j)
        D(i, j) =
            gf.kernelD(probe.normal().normalized(), source.center(), probe.center());
    }
  }

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(D, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }

  return D;
}

Eigen::MatrixXd anisotropicTEpsilon(const Cavity & cav, const IGreensFunction & gf_i,
                                    const IGreensFunction & gf_o) {
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  TIMER_ON("Computing SI");
  Eigen::MatrixXd SI = computeS(cav, gf_i);
  TIMER_OFF("Computing SI");
  TIMER_ON("Computing DI");
  Eigen::MatrixXd DI = computeD(cav, gf_i);
  TIMER_OFF("Computing DI");
  TIMER_ON("Computing SE");
  Eigen::MatrixXd SE = computeS(cav, gf_o);
  TIMER_OFF("Computing SE");
  TIMER_ON("Computing DE");
  Eigen::MatrixXd DE = computeD(cav, gf_o);
  TIMER_OFF("Computing DE");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  TIMER_ON("Assemble T matrix");
  Eigen::MatrixXd T = ((2 * M_PI * Id - DE * a) * SI +
                       SE * (2 * M_PI * Id + a * DI.adjoint().eval()));
  TIMER_OFF("Assemble T matrix");

  return T;
}

Eigen::MatrixXd isotropicTEpsilon(const Cavity & cav, const IGreensFunction & gf_i,
                                  double epsilon) {
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();

  // Compute SI and DI on the whole cavity, regardless of symmetry
  TIMER_ON("Computing SI");
  Eigen::MatrixXd SI = computeS(cav, gf_i);
  TIMER_OFF("Computing SI");
  TIMER_ON("Computing DI");
  Eigen::MatrixXd DI = computeD(cav, gf_i);
  TIMER_OFF("Computing DI");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  double fact = (epsilon + 1.0) / (epsilon - 1.0);
  TIMER_ON("Assemble T matrix");
  Eigen::MatrixXd T = (2 * M_PI * fact * Id - DI * a) * SI;
  TIMER_OFF("Assemble T matrix");

  return T;
}

Eigen::MatrixXd anisotropicRinfinity(const Cavity & cav,
                                     const IGreensFunction & gf_i,
                                     const IGreensFunction & gf_o) {
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  TIMER_ON("Computing SI");
  Eigen::MatrixXd SI = computeS(cav, gf_i);
  TIMER_OFF("Computing SI");
  TIMER_ON("Computing DI");
  Eigen::MatrixXd DI = computeD(cav, gf_i);
  TIMER_OFF("Computing DI");
  TIMER_ON("Computing SE");
  Eigen::MatrixXd SE = computeS(cav, gf_o);
  TIMER_OFF("Computing SE");
  TIMER_ON("Computing DE");
  Eigen::MatrixXd DE = computeD(cav, gf_o);
  TIMER_OFF("Computing DE");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  TIMER_ON("Assemble R matrix");
  Eigen::MatrixXd R =
      ((2 * M_PI * Id - DE * a) - SE * SI.llt().solve((2 * M_PI * Id - DI * a)));
  TIMER_OFF("Assemble R matrix");

  return R;
}

Eigen::MatrixXd isotropicRinfinity(const Cavity & cav,
                                   const IGreensFunction & gf_i) {
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();

  // Compute DI on the whole cavity, regardless of symmetry
  TIMER_ON("Computing DI");
  Eigen::MatrixXd DI = computeD(cav, gf_i);
  TIMER_OFF("Computing DI");

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  TIMER_ON("Assemble R matrix");
  Eigen::MatrixXd R = (2 * M_PI * Id - DI * a);
  TIMER_OFF("Assemble R matrix");

  return R;
}

Eigen::MatrixXd anisotropicIEFMatrix(const Cavity & cav,
                                     const IGreensFunction & gf_i,
                                     const IGreensFunction & gf_o) {
  Eigen::MatrixXd T = anisotropicTEpsilon(cav, gf_i, gf_o);
  Eigen::MatrixXd R = anisotropicRinfinity(cav, gf_i, gf_o);

  TIMER_ON("Assemble T^-1R matrix");
  Eigen::MatrixXd fullPCMMatrix = T.partialPivLu().solve(R);
  TIMER_OFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}

Eigen::MatrixXd isotropicIEFMatrix(const Cavity & cav, const IGreensFunction & gf_i,
                                   double epsilon) {
  Eigen::MatrixXd T = isotropicTEpsilon(cav, gf_i, epsilon);
  Eigen::MatrixXd R = isotropicRinfinity(cav, gf_i);

  TIMER_ON("Assemble T^-1R matrix");
  Eigen::MatrixXd fullPCMMatrix = T.partialPivLu().solve(R);
  TIMER_OFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}
} // namespace solver
