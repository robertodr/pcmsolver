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

#ifndef TDSOLVERDATA_HPP
#define TDSOLVERDATA_HPP

#include "Config.hpp"

/*! @struct TDSolverData
 *  @brief Contains all data defined from user input in the solver section.
 */

struct TDSolverData
{
    /*! The way the derivatives of the Green's function are evaluated */
    int howDerivative;
    /*! Evaluation policy for the diagonal of the S and D operators */
    int howIntegrator;
    /*! Static permittivity */
    double epsilonStatic;
    /*! Dynamic permittivity */
    double epsilonDynamic;
    /*! Relaxatione time */
    double tau;
    /*! The correction factor to be use in a CPCM calculation */
    double correction;
    /*! IEF relaxation time */
    double tauIEF;
    /*! Whether to use Cholesky decomposition */
    bool cholesky;
    /*! Whether the structure was initialized with user input or not */
    bool empty;

    TDSolverData() { empty = true; }
    TDSolverData(int how_d, int how_i,
        double es, double ed, double t, double corr, double tau, bool chol) :
        howDerivative(how_d), howIntegrator(how_i),
        epsilonStatic(es), epsilonDynamic(ed), tau(t),
        correction(corr), tauIEF(tau), cholesky(chol) { empty = false; }
};

#endif // SOLVERDATA_HPP
