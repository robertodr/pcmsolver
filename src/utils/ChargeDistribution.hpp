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

#ifndef CHARGEDISTRIBUTION_HPP
#define CHARGEDISTRIBUTION_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file ChargeDistribution.hpp
 *  \class ChargeDistribution
 *  \brief Class representing a classical charge distribution
 *  \author Roberto Di Remigio
 *  \date 2016
 */
class ChargeDistribution
{
  private:
    /*! Monopoles */
    Eigen::VectorXd monopoles_;
    /*! Monopoles sites */
    Eigen::Matrix3Xd monopolesSites_;
  public:
    ChargeDistribution(const Eigen::VectorXd & chg, const Eigen::Matrix3Xd & pos);
    /*! Return monopoles */
    Eigen::VectorXd monopoles() const { return monopoles_; }
    /*! Return monopoles sites */
    Eigen::Matrix3Xd monopolesSites() const { return monopolesSites_; }
    /*! Return i-th monopole */
    double monopoles(int i) const { return monopoles_(i); }
    /*! Return i-th monopole site */
    Eigen::Vector3d monopolesSites(int i) const { return monopolesSites_.col(i); }
};

/*! \typedef GreensFunctionValue
 *  \brief functor handle to the calculation of the value of a Greens's function in a point
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> GreensFunctionValue;

/*! \brief Computes Newton potential in a set of points for given Green's function and classical charge distribution
 *  \param[in] gf   the Green's function
 *  \param[in] grid where to evaluate the Newton potential
 *  \param[in] dist classical charge distribution
 *  \return the Newton potential on the grid
 */
Eigen::VectorXd computeNewtonPotential(const GreensFunctionValue & gf, const Eigen::Matrix3Xd & grid, const ChargeDistribution & dist);

#endif // CHARGEDISTRIBUTION_HPP
