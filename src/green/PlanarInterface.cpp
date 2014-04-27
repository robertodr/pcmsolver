/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#include "PlanarInterface.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double PlanarInterface<T>::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    throw std::runtime_error("Green's function for a planar interface has not yet been implemented!");
}

template<typename T>
T PlanarInterface<T>::operator()(T * sp, T * pp) const
{
    throw std::runtime_error("Green's function for a planar interface has not yet been implemented!");
}

template <typename T>
std::ostream & PlanarInterface<T>::printObject(std::ostream & os)
{
    os << "Green's function type: planar interface" << std::endl;
    os << "Permittivity (layer 1) = " << eps1_ << std::endl;
    os << "Permittivity (layer 2) = " << eps2_ << std::endl;
    os << "Position               = " << pos_ << std::endl;
    os << "Width                  = " << width_;
    return os;
}

template class PlanarInterface<double>;
template class PlanarInterface<AD_directional>;
template class PlanarInterface<AD_gradient>;
template class PlanarInterface<AD_hessian>;
