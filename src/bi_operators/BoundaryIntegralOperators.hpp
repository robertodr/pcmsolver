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

#ifndef BOUNDARYINTEGRALOPERATORS_HPP
#define BOUNDARYINTEGRALOPERATORS_HPP

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class Element;
class IGreensFunction;

namespace integrator {

class BoundaryIntegralOperator {
public:
  Eigen::MatrixXd operator()(const Cavity & cav, const IGreensFunction & gf) const;

private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const = 0;
};

class CollocationS __final : public BoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};

class NumericalS __final : public BoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};

class CollocationD __final : public BoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};

class PurisimaD __final : public BoundaryIntegralOperator {
private:
  /*! Computes the matrix representation of the double layer operator by collocation
   *  using the Purisima sum rule to compute the diagonal elements.
   *  \param[in] cav discretized cavity
   *  \param[in] gf  a Green's function
   *
   *  The sum rule for the diagonal elements is:
   *  \f[
   *    D_{ii} = -\left(2\pi + \sum_{j\neq i}D_{ij}a_j \right)\frac{1}{a_i}
   *  \f]
   */
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};

class NumericalD __final : public BoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};
} // namespace integrator

#endif // BOUNDARYINTEGRALOPERATORS_HPP
