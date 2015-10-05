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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef TDSOLVERHELPERFUNCTIONS_HPP
#define TDSOLVERHELPERFUNCTIONS_HPP

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

namespace td_solver {
inline Eigen::VectorXd K(const Eigen::VectorXd & Lambda, double factor) {
    Eigen::VectorXd K = Eigen::VectorXd::Zero(Lambda.size());
    for (int i = 0; i < Lambda.size(); ++i) {
        K(i) = (2 * M_PI - Lambda(i)) / (2 * M_PI * factor - Lambda(i));
    }
    return K;
}

inline Eigen::VectorXd tau(const Eigen::VectorXd & Lambda, double e_d, double e_0, double tau_D) {
    Eigen::VectorXd tau = Eigen::VectorXd::Zero(Lambda.size());
    for (int i = 0; i < Lambda.size(); ++i) {
        double num = (2 * M_PI - Lambda(i)) * e_d + 2 * M_PI + Lambda(i);
        double denom = (2 * M_PI - Lambda(i)) * e_0 + 2 * M_PI + Lambda(i);
        tau(i) = tau_D * num / denom;
    }
    return tau;
}

inline Eigen::VectorXd tauInverse(const Eigen::VectorXd & Lambda, double e_d, double e_0, double tau_D) {
    Eigen::VectorXd tau_inv = Eigen::VectorXd::Zero(Lambda.size());
    for (int i = 0; i < Lambda.size(); ++i) {
        double num = (2 * M_PI - Lambda(i)) * e_d + 2 * M_PI + Lambda(i);
        double denom = (2 * M_PI - Lambda(i)) * e_0 + 2 * M_PI + Lambda(i);
        tau_inv(i) = denom / (tau_D * num);
    }
    return tau_inv;
}
} // namespace td_solver

#endif // TDSOLVERHELPERFUNCTIONS_HPP
