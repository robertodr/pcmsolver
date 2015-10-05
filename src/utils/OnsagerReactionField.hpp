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

#ifndef ONSAGERREACTIONFIELD_HPP
#define ONSAGERREACTIONFIELD_HPP

#include <cmath>

/*! \brief Onsager reaction field for a point-like dipole in a spherical cavity
 *  \param[in] radius radius of the cavity
 *  \param[in] eps_0 static solvent permittivity
 *  \param[in] eps_d dynamic solvent permittivity
 *  \param[in] tau solvent relaxation time
 *  \param[in] t time step
 *  \note In \cite Corni2014 there's a typo in Eq. (50)
 */
double reactionField(double radius, double eps_0, double eps_d, double tau, double t)
{
    double tauOns = tau * (2 * eps_d + 1) / (2 * eps_0 + 1);
    double tmp_a = (2 * eps_d - 2) / (std::pow(radius, 3) * (2 * eps_d + 1));
    double tmp_b = (6 * (eps_0 - eps_d)) / ((2 * eps_d - 2) * (2 * eps_0 + 1));
    double tmp_c = std::exp(-t/tauOns) - 1;

    return tmp_a * (1 - tmp_b * tmp_c);
}

/**@{ Analytic expressions for a spherical cavity as in \cite Corni2014 */
/*! \brief Analytic Lambda matrix for spherical cavity
 *  \param[in] l angular momentum
 *  \note Eq. (15) in \cite Corni2014
 */
double Lambda_lm(int l)
{
    return -(2 * M_PI) / (2 * l + 1);
}

/*! \brief Analytic K matrix for spherical cavity
 *  \param[in] eps solvent permittivity
 *  \param[in] l angular momentum
 *  \note Eq. (16) in \cite Corni2014. It is related to the field factors in the multipolar models.
 */
double K_lm(double eps, int l)
{
    return (eps - 1) / (eps + l /(l + 1));
}

/*! \brief Analytic form of the relaxation times matrix
 *  \param[in] eps_0 static solvent permittivity
 *  \param[in] eps_d dynamic solvent permittivity
 *  \param[in] tau Debye relaxation time
 *  \param[in] l angular momentum
 *  \note Eq. (33) in \cite Corni2014
 */
double tau_lm(double eps_0, double eps_d, double tau, int l)
{
    return tau * ((l + 1) * eps_d + l) / ((l + 1) * eps_0 + l);
}
/**@}*/

#endif // ONSAGERREACTIONFIELD_HPP
