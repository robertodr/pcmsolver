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

#include "VIEFSolver.hpp"

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

void VIEFSolver::buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  isotropic_ = (gf_i.uniform() && gf_o.uniform());
  isotropic_ ? buildIsotropicMatrix(cavity, gf_i, gf_o) : buildAnisotropicMatrix(cavity, gf_i, gf_o);
}

void VIEFSolver::buildAnisotropicMatrix(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  R_infinity_ = anisotropicRinfinity(cav, gf_i, gf_o);
  PCMMatrix_ = anisotropicTEpsilon(cav, gf_i, gf_o) * R_infinity_.adjoint().eval();
  // Symmetrize K := (K + K+)/2
  hermitivitize(PCMMatrix_);
  // Pack into a block diagonal matrix
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  // For the moment just packs into a std::vector<Eigen::MatrixXd>
  symmetryPacking(blockPCMMatrix_, PCMMatrix_, dimBlock, nrBlocks);
  symmetryPacking(blockR_infinity_, R_infinity_, dimBlock, nrBlocks);
  built_ = true;
}

void VIEFSolver::buildIsotropicMatrix(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  R_infinity_ = isotropicRinfinity(cav, gf_i);
  PCMMatrix_ = isotropicTEpsilon(cav, gf_i, profiles::epsilon(gf_o.permittivity())) * R_infinity_.adjoint().eval();
  // Symmetrize K := (K + K+)/2
  hermitivitize(PCMMatrix_);
  // Pack into a block diagonal matrix
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  // For the moment just packs into a std::vector<Eigen::MatrixXd>
  symmetryPacking(blockPCMMatrix_, PCMMatrix_, dimBlock, nrBlocks);
  symmetryPacking(blockR_infinity_, R_infinity_, dimBlock, nrBlocks);

  built_ = true;
}

inline Eigen::VectorXd VIEFSolver::getBareASC(const Eigen::VectorXd & dressedASC, int irrep) const
{
  int irrDim = blockR_infinity_[irrep].rows();
  Eigen::VectorXd bareASC = Eigen::VectorXd::Zero(dressedASC.size());
  bareASC.segment(irrep*irrDim, irrDim) = blockR_infinity_[irrep].adjoint() * dressedASC.segment(irrep*irrDim, irrDim);
  return bareASC;
}

inline Eigen::VectorXd VIEFSolver::getDressedMEP(const Eigen::VectorXd & bareMEP, int irrep) const
{
  int irrDim = blockR_infinity_[irrep].rows();
  Eigen::VectorXd dressedMEP = Eigen::VectorXd::Zero(bareMEP.size());
  dressedMEP.segment(irrep*irrDim, irrDim) = blockR_infinity_[irrep] * bareMEP.segment(irrep*irrDim, irrDim);
  return dressedMEP;
}

Eigen::VectorXd VIEFSolver::computeCharge_impl(const Eigen::VectorXd & potential, int irrep, double CGtol) const
{
  // The potential and charge vector are of dimension equal to the
  // full dimension of the cavity. We have to select just the part
  // relative to the irrep needed.
  int fullDim = PCMMatrix_.rows();
  Eigen::VectorXd ASC = Eigen::VectorXd::Zero(fullDim);
  int nrBlocks = blockPCMMatrix_.size();
  int irrDim = fullDim/nrBlocks;
  // Initialize Conjugate Gradient solver
  // use default maximum number of iterations and tolerance
  Eigen::ConjugateGradient<Eigen::MatrixXd> CGSolver;
  CGSolver.compute(blockPCMMatrix_[irrep]);
  CGSolver.setTolerance(CGtol);
  // Preprocess incoming potential, get only the relevant irrep
  Eigen::VectorXd tildeMEP = getDressedMEP(potential, irrep);
  // Obtain \tilde{q} by solving \tilde{Y}\tilde{q} + \tilde{v} = 0 only for the relevant irrep
  ASC.segment(irrep*irrDim, irrDim) = CGSolver.solve(-tildeMEP.segment(irrep*irrDim, irrDim));
  return getBareASC(ASC, irrep);
}

Eigen::VectorXd VIEFSolver::initialGuess_impl(const Eigen::VectorXd & MEP, double nuc_chg, int irrep) const
{
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(MEP.size());
  switch(guess_) {
    case Trivial:     /* Do nothing */
                      break;
    case Uniform:     tmp = initialGuessUniform(nuc_chg, irrep);
                      break;
    case Diagonal:    tmp = getBareASC(initialGuessDiagonal(getDressedMEP(MEP, irrep), irrep), irrep);
                      break;
    case LowAccuracy: tmp = initialGuessLowAccuracy(MEP, irrep);
                      break;
  }
  return tmp;
}

Eigen::VectorXd VIEFSolver::error_impl(const Eigen::VectorXd & dressedASC,
    const Eigen::VectorXd & bareMEP, int irrep) const
{
  int fullDim = PCMMatrix_.rows();
  int nrBlocks = blockPCMMatrix_.size();
  int irrDim = fullDim/nrBlocks;
  Eigen::VectorXd error = Eigen::VectorXd::Zero(fullDim);
  error.segment(irrep*irrDim, irrDim) = blockPCMMatrix_[irrep] * dressedASC.segment(irrep*irrDim, irrDim)
    + getDressedMEP(bareMEP, irrep).segment(irrep*irrDim, irrDim);
  return error;
}

Eigen::VectorXd VIEFSolver::updateCharge_impl(const Eigen::VectorXd & dressedASC,
    const Eigen::VectorXd & residual, int irrep) const
{
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(dressedASC.size());
  switch(update_) {
    case SSD:        tmp = updateChargeSSD(dressedASC, residual, irrep);
                     break;
    case LineSearch: tmp = updateChargeLineSearch(dressedASC, residual, irrep);
                     break;
  }
  return getBareASC(tmp, irrep);
}

std::ostream & VIEFSolver::printSolver(std::ostream & os)
{
  std::string type;
  if (isotropic_) {
    type = "Variational IEFPCM, isotropic";
  } else {
    type = "Variational IEFPCM, anisotropic";
  }
  os << "Solver Type:       " << type << std::endl;
  os << "ASC initial guess: " << guess(guess_) << std::endl;
  os << "ASC update:        " << update(update_);
  return os;
}

