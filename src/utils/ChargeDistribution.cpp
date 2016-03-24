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

#include "ChargeDistribution.hpp"

#include "Config.hpp"

#include <Eigen/Core>

ChargeDistribution::ChargeDistribution(const Eigen::VectorXd & chg,
    const Eigen::Matrix3Xd & pos) : monopoles_(chg), monopolesSites_(pos) {};

Eigen::VectorXd computeNewtonPotential(const GreensFunctionValue & gf, const Eigen::Matrix3Xd & grid, const ChargeDistribution & dist)
{
  Eigen::VectorXd newton = Eigen::VectorXd::Zero(grid.cols());
  int nSites = dist.monopoles().size();
  for (int i = 0; i < nSites; ++i) {
    for (int j = 0; j < grid.cols(); ++j) {
      newton(j) += dist.monopoles(i) * gf(grid.col(j), dist.monopolesSites(i));
    }
  }
  return newton;
}
