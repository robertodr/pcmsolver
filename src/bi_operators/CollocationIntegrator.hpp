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

#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Element;
class IGreensFunction;

#include "BoundaryIntegralOperator.hpp"

/*! \file CollocationIntegrator.hpp
 *  \struct CollocationIntegrator
 *  \brief Implementation of the single and double layer operators matrix
 *representation using one-point collocation
 *  \author Roberto Di Remigio
 *  \date 2015, 2016
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *   S_{ii} = factor * \sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *   D_{ii} = -factor * \sqrt{\frac{\pi}{a_i}} \frac{1}{R_I}
 *  \f]
 */

namespace integrator {
class CollocationS __final : public BoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};

class CollocationD __final : public BoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd compute(const std::vector<Element> & elems,
                                  const IGreensFunction & gf) const __override;
};
} // namespace integrator

#endif // COLLOCATIONINTEGRATOR_HPP
