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

#include "TDIEFSolver.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "bi_operators/IBoundaryIntegralOperator.hpp"
#include "cavity/ICavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "TDSolverData.hpp"
#include "Debye.hpp"
#include "ITDSolver.hpp"

namespace pcm {
namespace td_solver {
TDIEFSolver::TDIEFSolver(double es, double ed, double t, bool cholesky)
    : ITDSolver(es, ed, t), useCholesky_(cholesky) {}

void TDIEFSolver::buildSystemMatrix_impl(const ICavity & cavity,
                                         const IGreensFunction & gf_i,
                                         const IBoundaryIntegralOperator & op) {
  useCholesky_ ? systemMatrix_Cholesky(cavity, gf_i, op)
               : systemMatrix_Lowdin(cavity, gf_i, op);
}

void TDIEFSolver::systemMatrix_Lowdin(const ICavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const IBoundaryIntegralOperator & op) {
  using namespace td_solver;
  // The total size of the cavity
  int cavitySize = cavity.size();
  // Diagonal areas matrix
  Eigen::MatrixXd A = cavity.elementArea().asDiagonal();
  // Identiy matrix
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
  double e_d = permittivity_.epsilonDynamic;
  double e_0 = permittivity_.epsilonStatic;

  // Compute S on the whole cavity, regardless of symmetry
  Eigen::MatrixXd S = op.computeS(cavity, gf_i);
  // Compute D on whole cavity, regardless of symmetry
  Eigen::MatrixXd D = op.computeD(cavity, gf_i);
  // Make sure S is symmetric (might be unnecessary)
  utils::hermitivitize(S);
  // Compute S^1/2 and S^-1/2
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> S_eigen(S);
  Eigen::MatrixXd S_sqrt = S_eigen.operatorSqrt();
  Eigen::MatrixXd S_invsqrt = S_eigen.operatorInverseSqrt();
  // Compute S^-1/2DAS^1/2 and symmetrize it
  D = S_invsqrt * D * A * S_sqrt;
  utils::hermitivitize(D);
  // Diagonalize S^-1/2DAS^1/2
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> D_eigen(D);
  // Get Lambda and the transformation matrix
  Lambda_ = D_eigen.eigenvalues();
  Eigen::MatrixXd T = D_eigen.eigenvectors();
  // Form tau^-1 (relaxation times matrix)
  tau_ = detail::tau(Lambda_, e_d, e_0, permittivity_.tau);
  Eigen::VectorXd tau_inv = detail::tauInverse(Lambda_, e_d, e_0, permittivity_.tau);
  // Form A_ (the dynamic matrix)
  double f_d = (e_d + 1.0) / (e_d - 1.0);
  K_d_ = detail::K(Lambda_, f_d);
  A_ = -S_invsqrt * T * K_d_.asDiagonal() * T.adjoint().eval() * S_invsqrt;
  // Form B_ (the stati matrix)
  double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
  K_0_ = detail::K(Lambda_, f_0);
  B_ = -S_invsqrt * T * tau_inv.asDiagonal() * K_0_.asDiagonal() *
       T.adjoint().eval() * S_invsqrt;
  // Form C_ (the relaxation times matrix)
  C_ = S_invsqrt * T * tau_inv.asDiagonal() * T.adjoint().eval() * S_sqrt;
  // Form PCMMatrix_ (needed to form the initial value of the ASC)
  PCMMatrix_ = -S_invsqrt * T * K_d_.asDiagonal() * T.adjoint().eval() * S_invsqrt;

  built_ = true;
}

void TDIEFSolver::systemMatrix_Cholesky(const ICavity & /* cavity */,
                                        const IGreensFunction & /* gf_i */,
                                        const IBoundaryIntegralOperator & /* op */) {
  /*
  using namespace td_solver;
  // The total size of the cavity
  int cavitySize = cavity.size();
  // Diagonal areas matrix
  Eigen::MatrixXd A = cavity.elementArea().asDiagonal();
  // Identiy matrix
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
  double e_d = permittivity_.epsilonDynamic;
  double e_0 = permittivity_.epsilonStatic;

  Vacuum<DerivativeTraits, IntegratorPolicy> gf_i;
  // Compute S on the whole cavity, regardless of symmetry
  Eigen::MatrixXd S = gf_i.singleLayer(cavity.elements());
  // Compute D on whole cavity, regardless of symmetry
  Eigen::MatrixXd D = gf_i.doubleLayer(cavity.elements());
  // Make sure S is symmetric (might be unnecessary)
  utils::hermitivitize(S);
  // Compute Cholesky decomposition of S = LVL^t
  Eigen::LDLT<Eigen::MatrixXd> S_lvlt(S);
  Eigen::MatrixXd L = S_lvlt.matrixL();
  Eigen::MatrixXd L_inv = L.inverse();
  Eigen::MatrixXd Lt_inv = L_inv.transpose();
  Eigen::MatrixXd V = S_lvlt.vectorD().asDiagonal();
  Eigen::MatrixXd V_inv = V.inverse();
  // Compute L^-1DAL and symmetrize it
  D = L.inverse() * D * A * L;
  //utils::hermitivitize(D);
  // Diagonalize L^-1DAL
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> D_eigen(D);
  // Get Lambda and the transformation matrix
  Eigen::MatrixXd Lambda = D_eigen.eigenvalues().asDiagonal();
  std::cout << D_eigen.eigenvalues() << std::endl;
  Eigen::MatrixXd T = D_eigen.eigenvectors();
  // Form tau^-1 (relaxation times matrix)
  Eigen::MatrixXd tau_inv = detail::tauInverse(Lambda, e_d, e_0, permittivity_.tau);
  // Form A_ (the dynamic matrix)
  double f_d = (e_d + 1.0) / (e_d - 1.0);
  A_ = - Lt_inv * T * detail::K(Lambda, f_d) * T.adjoint().eval() * V_inv * L_inv;
  // Form B_ (the static matrix)
  double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
  B_ = - Lt_inv * T * tau_inv * detail::K(Lambda, f_0) * T.adjoint().eval() * V_inv *
  L_inv;
  // Form C_ (the relaxation times matrix)
  C_ = Lt_inv * T * tau_inv * T.adjoint().eval() * L * V;
  // Form PCMMatrix_ (needed to form the initial value of the ASC)
  PCMMatrix_ = - Lt_inv * T * detail::K(Lambda, f_d) * T.adjoint().eval() * V_inv *
  L_inv;

  built_ = true;
  */
}

Eigen::VectorXd TDIEFSolver::propagateASC_impl(
    double dt,
    const Eigen::VectorXd & MEP_current,
    const Eigen::VectorXd & MEP_previous,
    const Eigen::VectorXd & ASC_previous) const {
  return (A_ * (MEP_current - MEP_previous) + dt * B_ * MEP_previous -
          dt * C_ * ASC_previous + ASC_previous);
}

Eigen::VectorXd TDIEFSolver::initialValueASC_impl(
    const Eigen::VectorXd & MEP) const {
  return (PCMMatrix_ * MEP);
}

std::ostream & TDIEFSolver::printSolver(std::ostream & os) {
  os << "Solver Type: IEFPCM, isotropic" << std::endl;
  os << permittivity_;
  return os;
}

ITDSolver * createTDIEFSolver(const TDSolverData & data) {
  return new TDIEFSolver(
      data.epsilonStatic, data.epsilonDynamic, data.tau, data.cholesky);
}

namespace detail {
Eigen::VectorXd K(const Eigen::VectorXd & Lambda, double factor) {
  Eigen::VectorXd K = Eigen::VectorXd::Zero(Lambda.size());
  for (int i = 0; i < Lambda.size(); ++i) {
    K(i) = (2 * M_PI - Lambda(i)) / (2 * M_PI * factor - Lambda(i));
  }
  return K;
}

Eigen::VectorXd tau(const Eigen::VectorXd & Lambda,
                    double e_d,
                    double e_0,
                    double tau_D) {
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(Lambda.size());
  for (int i = 0; i < Lambda.size(); ++i) {
    double num = (2 * M_PI - Lambda(i)) * e_d + 2 * M_PI + Lambda(i);
    double denom = (2 * M_PI - Lambda(i)) * e_0 + 2 * M_PI + Lambda(i);
    tau(i) = tau_D * num / denom;
  }
  return tau;
}

Eigen::VectorXd tauInverse(const Eigen::VectorXd & Lambda,
                           double e_d,
                           double e_0,
                           double tau_D) {
  Eigen::VectorXd tau_inv = Eigen::VectorXd::Zero(Lambda.size());
  for (int i = 0; i < Lambda.size(); ++i) {
    double num = (2 * M_PI - Lambda(i)) * e_d + 2 * M_PI + Lambda(i);
    double denom = (2 * M_PI - Lambda(i)) * e_0 + 2 * M_PI + Lambda(i);
    tau_inv(i) = denom / (tau_D * num);
  }
  return tau_inv;
}
} // namespace detail
} // namespace td_solver
} // namespace pcm
