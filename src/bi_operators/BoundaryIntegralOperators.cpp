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

#include "BoundaryIntegralOperators.hpp"

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/map.hpp>

#include "IntegratorHelperFunctions.hpp"
#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "utils/QuadratureRules.hpp"

namespace integrator {
Eigen::MatrixXd BoundaryIntegralOperator::operator()(
    const Cavity & cav, const IGreensFunction & gf) const {
  Eigen::MatrixXd biop = compute(cav.elements(), gf);
  // Perform symmetry blocking
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(biop, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }
  return biop;
}

Eigen::MatrixXd CollocationS::compute(const std::vector<Element> & elems,
                                      const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    S(i, i) = gf.singleLayer(source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        S(i, j) = gf.kernelS(source.center(), probe.center());
    }
  }
  return S;
}

Eigen::MatrixXd NumericalS::compute(const std::vector<Element> & elems,
                                    const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    S(i, i) = integrateS<32, 16>(gf.exportKernelS(), source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        S(i, j) = gf.kernelS(source.center(), probe.center());
    }
  }
  return S;
}

Eigen::MatrixXd CollocationD::compute(const std::vector<Element> & elems,
                                      const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    D(i, i) = gf.doubleLayer(source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        D(i, j) =
            gf.kernelD(probe.normal().normalized(), source.center(), probe.center());
    }
  }
  return D;
}

Eigen::MatrixXd PurisimaD::compute(const std::vector<Element> & elems,
                                   const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    double D_ii = 0.0;
    Element source = elems[i];
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j) {
        D(i, j) =
            gf.kernelD(probe.normal().normalized(), source.center(), probe.center());
        D_ii += D(i, j) * probe.area();
      }
    }
    D(i, i) = -(2 * M_PI + D_ii) / (source.area());
  }
  return D;
}

Eigen::MatrixXd NumericalD::compute(const std::vector<Element> & elems,
                                    const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    D(i, i) = integrateD<32, 16>(gf.exportKernelD(), source);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        D(i, j) =
            gf.kernelD(probe.normal().normalized(), source.center(), probe.center());
    }
  }
  return D;
}
} // namespace integrator
