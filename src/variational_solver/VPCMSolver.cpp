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

#include "VPCMSolver.hpp"

#include <iostream>
#include <cmath>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "MathUtils.hpp"

Eigen::VectorXd VPCMSolver::updateChargeSSD(const Eigen::VectorXd & dressedASC,
    const Eigen::VectorXd & bareMEP, int irrep) const
{
  int fullDim = PCMMatrix_.rows();
  int nrBlocks = blockPCMMatrix_.size();
  int irrDim = fullDim/nrBlocks;
  // Get the residual
  Eigen::VectorXd residual = -error_impl(dressedASC, bareMEP, irrep);

  Eigen::VectorXd updated = Eigen::VectorXd::Zero(fullDim);
  updated.segment(irrep*irrDim, irrDim) = dressedASC.segment(irrep*irrDim, irrDim)
    + (residual.segment(irrep*irrDim, irrDim)).cwiseQuotient(blockPCMMatrix_[irrep].diagonal());
  return updated;
}

Eigen::VectorXd VPCMSolver::updateChargeLineSearch(const Eigen::VectorXd & dressedASC,
    const Eigen::VectorXd & bareMEP, int irrep) const
{
  int fullDim = PCMMatrix_.rows();
  int nrBlocks = blockPCMMatrix_.size();
  int irrDim = fullDim/nrBlocks;
  // Get the residual
  Eigen::VectorXd residual = -error_impl(dressedASC, bareMEP, irrep);
  // Calculate the coefficient:
  double num = (residual.segment(irrep*irrDim, irrDim)).squaredNorm();
  double den = (residual.segment(irrep*irrDim, irrDim).transpose()
      * blockPCMMatrix_[irrep] * residual.segment(irrep*irrDim, irrDim));
  double alpha = num / den;

  Eigen::VectorXd updated = Eigen::VectorXd::Zero(fullDim);
  updated.segment(irrep*irrDim, irrDim) =
    dressedASC.segment(irrep*irrDim, irrDim) + alpha * residual.segment(irrep*irrDim, irrDim);
  return updated;
}
