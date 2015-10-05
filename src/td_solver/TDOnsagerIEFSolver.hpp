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

#ifndef TDONSAGERIEFSOLVER_HPP
#define TDONSAGERIEFSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "Cavity.hpp"
#include "Debye.hpp"
#include "Element.hpp"
#include "MathUtils.hpp"
#include "Vacuum.hpp"
#include "TDSingleIEFSolver.hpp"

/*! \file TDOnsagerIEFSolver.hpp
 *  \class TDOnsagerIEFSolver
 *  \brief Time-dependent solver for isotropic IEF with a single Onsager relaxation time
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits,
          typename IntegratorPolicy>
class TDOnsagerIEFSolver : public TDSingleIEFSolver<DerivativeTraits, IntegratorPolicy>
{
public:
    TDOnsagerIEFSolver() {}
    /*! \brief Construct solver
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     */
    TDOnsagerIEFSolver(double es, double ed, double t)
            : TDSingleIEFSolver<DerivativeTraits, IntegratorPolicy>(es, ed, t, (t*(2*ed+1)/(2*es+1))) {}
    virtual ~TDOnsagerIEFSolver() {}
    friend std::ostream & operator<<(std::ostream & os, TDOnsagerIEFSolver & solver) {
        return solver.printSolver(os);
    }
private:
    virtual std::ostream & printSolver(std::ostream & os) {
        os << "Solver Type: IEFPCM, isotropic" << std::endl;
        os << "Onsager relaxation time = " << this->tauIEF_;
        return os;
    }
};

#endif // TDONSAGERIEFSOLVER_HPP
