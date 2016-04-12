/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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

#ifndef DEBYE_HPP
#define DEBYE_HPP

#include <cmath>
#include <iosfwd>

#include "Config.hpp"

/*! \file Debye.hpp
 *  \struct Debye
 *  \brief A time-dependent Debye-type dielectric profile
 *  \author Roberto Di Remigio
 *  \date 2015
 */

struct Debye
{
    Debye() {}
    Debye(double es, double ed, double t) : epsilonStatic(es), epsilonDynamic(ed), tau(t) {}
    /// Static dielectric constant
    double epsilonStatic;
    /// Dynamic dielectric constant
    double epsilonDynamic;
    /// Relaxation time
    double tau;
    /*! Return static Onsager factor
     *  \param[in] corr correction factor
     */
    double staticOnsager(double corr = 0.0) const { return -(epsilonStatic - 1) / (epsilonStatic + corr); }
    /*! Return dynamic Onsager factor
     *  \param[in] corr correction factor
     */
    double dynamicOnsager(double corr = 0.0) const { return -(epsilonDynamic - 1) / (epsilonDynamic + corr); }
    friend std::ostream & operator<<(std::ostream & os, Debye & obj)
    {
        os << "Profile functional form: Debye" << std::endl;
        os << "Static permittivity  = " << obj.epsilonStatic << std::endl;
        os << "Dynamic permittivity = " << obj.epsilonDynamic << std::endl;
        os << "Relaxation time      = " << obj.tau * AUToFemtoseconds() << " fs";
        return os;
    }
};

#endif // DEBYE_HPP
