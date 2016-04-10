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

#ifndef TDPCMSOLVER_HPP
#define TDPCMSOLVER_HPP

#include <iosfwd>
#include <sstream>
#include <string>

#include "Config.hpp"

class Cavity;

#include "Debye.hpp"
#include "ErrorHandling.hpp"

/*! \file TDPCMSolver.hpp
 *  \class TDPCMSolver
 *  \brief Abstract Base Class for time-dependent solvers inheritance hierarchy.
 *  \author Roberto Di Remigio
 *  \date 2015
 *  We use the Non-Virtual Interface idiom.
 */

class TDPCMSolver
{
public:
    TDPCMSolver() {}
    /*! \brief Construct solver from two Green's functions
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     */
    TDPCMSolver(double es, double ed, double t)
            : permittivity_(Debye(es, ed, t)), built_(false) {}
    virtual ~TDPCMSolver() {}

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    void buildSystemMatrix(const Cavity & cavity) { buildSystemMatrix_impl(cavity); }
    /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
     *  \param[in] dt the time step for the Euler integrator
     *  \param[in] MEP_current the vector containing the MEP at cavity points, at time (t + dt)
     *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time t
     *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time t
     */
    Eigen::VectorXd propagateASC(double dt, const Eigen::VectorXd & MEP_current, const Eigen::VectorXd & MEP_previous,
            const Eigen::VectorXd & ASC_previous) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
        return propagateASC_impl(dt, MEP_current, MEP_previous, ASC_previous);
    }
    /*! \brief Returns the ASC at initial time
     *  \param[in] MEP the vector containing the MEP at cavity points at initial time
     *  Uses dynamic PCM matrix to compute the initial value for the ASC
     */
    Eigen::VectorXd initialValueASC(const Eigen::VectorXd & MEP) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
        return initialValueASC_impl(MEP);
    }
    std::string printEnvironment() {
        std::stringstream tmp;
        tmp << ".... Inside " << std::endl;
        tmp << "Green's function type: vacuum" << std::endl;
        tmp << ".... Outside " << std::endl;
        tmp << "Green's function type: uniform dielectric" << std::endl;
        tmp << permittivity_;
        return tmp.str();
    }

    friend std::ostream & operator<<(std::ostream & os, TDPCMSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    /*! Time-dependent isotropic, uniform Debye permittivity profile */
    Debye permittivity_;
    /*! Whether the system matrices have been built */
    bool built_;
    /*! Dynamic matrix */
    Eigen::MatrixXd A_;
    /*! Static matrix scaled by relaxation times */
    Eigen::MatrixXd B_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity) = 0;
    /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
     *  \param[in] dt the time step for the Euler integrator
     *  \param[in] MEP_current the vector containing the MEP at cavity points, at time (t + dt)
     *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time t
     *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time t
     */
    virtual Eigen::VectorXd propagateASC_impl(double dt, const Eigen::VectorXd & MEP_current, const Eigen::VectorXd & MEP_previous,
            const Eigen::VectorXd & ASC_previous) const = 0;
    /*! \brief Returns the ASC at initial time
     *  \param[in] MEP the vector containing the MEP at cavity points at initial time
     *  Uses dynamic PCM matrix to compute the initial value for the ASC
     */
    virtual Eigen::VectorXd initialValueASC_impl(const Eigen::VectorXd & MEP) const = 0;
    virtual std::ostream & printSolver(std::ostream & os) = 0;
};

#endif // TDPCMSOLVER_HPP
