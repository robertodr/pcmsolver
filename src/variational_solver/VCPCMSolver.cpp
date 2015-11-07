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

#include "VCPCMSolver.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>

#include "Cavity.hpp"
#include "Element.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"
#include "VSolverImpl.hpp"

void VCPCMSolver::buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  double epsilon = profiles::epsilon(gf_o.permittivity());
  // This is the inverse of the factor used in the
  // traditional CPCM solver!
  double fact = (epsilon + correction_)/(epsilon - 1);
  S_ = fact * gf_i.singleLayer(cavity.elements());
  hermitivitize(S_);

  // Symmetry-pack
  // The number of irreps in the group
  int nrBlocks = cavity.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cavity.irreducible_size();
  symmetryPacking(blockS_, S_, dimBlock, nrBlocks);

  built_ = true;
}

Eigen::VectorXd VCPCMSolver::computeCharge_impl(const Eigen::VectorXd & potential, int irrep) const
{
  // The potential and charge vector are of dimension equal to the
  // full dimension of the cavity. We have to select just the part
  // relative to the irrep needed.
  int fullDim = S_.rows();
  Eigen::VectorXd ASC = Eigen::VectorXd::Zero(fullDim);
  int nrBlocks = blockS_.size();
  int irrDim = fullDim/nrBlocks;
  // Initialize Conjugate Gradient solver
  // use default maximum number of iterations and tolerance
  Eigen::ConjugateGradient<Eigen::MatrixXd> CGSolver;
  CGSolver.compute(blockS_[irrep]);
  // Obtain q by solving \frac{1}{f(\varepsilon)}Sq + v = 0 only for the relevant irrep
  ASC.segment(irrep*irrDim, irrDim) = CGSolver.solve(-potential.segment(irrep*irrDim, irrDim));
  return ASC;
}

Eigen::VectorXd VCPCMSolver::updateCharge_impl(const Eigen::VectorXd & potential, int irrep) const
{}

std::ostream & VCPCMSolver::printSolver(std::ostream & os)
{
  os << "Solver Type: Variational C-PCM" << std::endl;
  return os;
}

