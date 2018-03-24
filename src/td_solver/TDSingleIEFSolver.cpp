/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

#include "TDSingleIEFSolver.hpp"

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "Debye.hpp"
#include "ITDSolver.hpp"
#include "TDSolverData.hpp"
#include "bi_operators/IBoundaryIntegralOperator.hpp"
#include "cavity/ICavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"

namespace pcm {
namespace td_solver {
TDSingleIEFSolver::TDSingleIEFSolver(double es, double ed, double t, double tau)
    : ITDSolver(es, ed, t), tauIEF_(tau) {}

void TDSingleIEFSolver::buildSystemMatrix_impl(
    const ICavity & cavity,
    const IGreensFunction & gf_i,
    const IBoundaryIntegralOperator & op) {
  // The total size of the cavity
  int cavitySize = cavity.size();
  // Identity matrix
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
  double e_d = permittivity_.epsilonDynamic;
  double e_0 = permittivity_.epsilonStatic;
  // Compute S and D on the whole cavity, regardless of symmetry
  Eigen::MatrixXd S = op.computeS(cavity, gf_i);
  Eigen::MatrixXd D = op.computeD(cavity, gf_i);
  Eigen::MatrixXd A = cavity.elementArea().asDiagonal();
  // Form A_ (the dynamic matrix)
  double f_d = (e_d + 1.0) / (e_d - 1.0);
  A_ = 2 * M_PI * f_d * S - D * A * S;
  Eigen::FullPivLU<Eigen::MatrixXd> A__LU(A_);
  if (!(A__LU.isInvertible()))
    PCMSOLVER_ERROR("A_ matrix is not invertible!");
  A_ = -A__LU.inverse();
  A_ *= (2 * M_PI * Id - D * A);
  utils::hermitivitize(A_);
  // Form B_ (the static matrix)
  double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
  B_ = 2 * M_PI * f_0 * S - D * A * S;
  Eigen::FullPivLU<Eigen::MatrixXd> B__LU(B_);
  if (!(B__LU.isInvertible()))
    PCMSOLVER_ERROR("B_ matrix is not invertible!");
  B_ = -B__LU.inverse();
  B_ *= (2 * M_PI * Id - D * A);
  utils::hermitivitize(B_);
  built_ = true;
}

Eigen::VectorXd TDSingleIEFSolver::propagateASC_impl(
    double dt,
    const Eigen::VectorXd & MEP_current,
    const Eigen::VectorXd & MEP_previous,
    const Eigen::VectorXd & ASC_previous) const {
  double factor = dt / tauIEF_;
  return (A_ * (MEP_current - MEP_previous) + factor * B_ * MEP_previous -
          factor * ASC_previous + ASC_previous);
}

Eigen::VectorXd TDSingleIEFSolver::initialValueASC_impl(
    const Eigen::VectorXd & MEP) const {
  return (A_ * MEP);
}

std::ostream & TDSingleIEFSolver::printSolver(std::ostream & os) {
  os << "Solver Type: IEFPCM, isotropic" << std::endl;
  os << "IEF relaxation time = " << tauIEF_;
  return os;
}

ITDSolver * createTDSingleIEFSolver(const TDSolverData & data) {
  return new TDSingleIEFSolver(
      data.epsilonStatic, data.epsilonDynamic, data.tau, data.tauIEF);
}

ITDSolver * createTDOnsagerIEFSolver(const TDSolverData & data) {
  double tauOnsager =
      data.tau * (2 * data.epsilonDynamic + 1) / (2 * data.epsilonStatic + 1);
  return new TDSingleIEFSolver(
      data.epsilonStatic, data.epsilonDynamic, data.tau, tauOnsager);
}
} // namespace td_solver
} // namespace pcm
