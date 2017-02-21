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

#include "TDCPCMSolver.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "Config.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "bi_operators/BoundaryIntegralOperator.hpp"
#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "TDSolverData.hpp"
#include "utils/Factory.hpp"
#include "Debye.hpp"

TDCPCMSolver::TDCPCMSolver(double es, double ed, double t, double corr)
    : TDPCMSolver(es, ed, t), correction_(corr) {}

void TDCPCMSolver::buildSystemMatrix_impl(const Cavity & cavity,
                                          const IGreensFunction & gf_i,
                                          const BoundaryIntegralOperator & op) {
  // Compute SI on the whole cavity, regardless of symmetry
  double f_d = permittivity_.dynamicOnsager(correction_);
  A_ = f_d * op.computeS(cavity, gf_i).partialPivLu().inverse();
  hermitivitize(A_);

  double f_0 = permittivity_.staticOnsager(correction_);
  B_ = f_0 * op.computeS(cavity, gf_i).partialPivLu().inverse();
  hermitivitize(B_);

  built_ = true;
}

Eigen::VectorXd TDCPCMSolver::propagateASC_impl(
    double dt,
    const Eigen::VectorXd & MEP_current,
    const Eigen::VectorXd & MEP_previous,
    const Eigen::VectorXd & ASC_previous) const {
  double tau_CPCM = permittivity_.tau *
                    (permittivity_.epsilonDynamic / permittivity_.epsilonStatic);
  double factor = dt / tau_CPCM;
  return (A_ * (MEP_current - MEP_previous) + factor * B_ * MEP_previous -
          factor * ASC_previous + ASC_previous);
}

Eigen::VectorXd TDCPCMSolver::initialValueASC_impl(
    const Eigen::VectorXd & MEP) const {
  return A_ * MEP;
}

std::ostream & TDCPCMSolver::printSolver(std::ostream & os) {
  os << "Solver Type: C-PCM" << std::endl;
  os << "Correction factor = " << correction_;
  return os;
}

namespace {
TDPCMSolver * createTDCPCMSolver(const TDSolverData & data) {
  return new TDCPCMSolver(
      data.epsilonStatic, data.epsilonDynamic, data.tau, data.correction);
}
const std::string TDCPCMSOLVER("TDCPCM");
const bool registeredTDCPCMSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDCPCMSOLVER,
        createTDCPCMSolver);
}
