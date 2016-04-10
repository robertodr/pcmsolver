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

#ifndef TDPCMITERATIVESOLVER_HPP
#define TDPCMITERATIVESOLVER_HPP

#include <iosfwd>
#include <sstream>
#include <string>

#include "Config.hpp"

class Cavity;

#include "Debye.hpp"
#include "ErrorHandling.hpp"

/*! \file TDPCMIterativeSolver.hpp
 *  \class TDPCMIterativeSolver
 *  \brief Abstract Base Class for time-dependent solvers inheritance hierarchy.
 *  \author Roberto Di Remigio
 *  \date 2015
 *  We use the Non-Virtual Interface idiom.
 */

class TDPCMIterativeSolver
{
public:
    TDPCMIterativeSolver() {}
    /*! \brief Construct solver
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     */
    TDPCMIterativeSolver(double es, double ed, double t)
            : permittivity_(Debye(es, ed, t)), built_(false), cgBuilt_(false) {}
    virtual ~TDPCMIterativeSolver() {}

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    void buildSystemMatrix(const Cavity & cavity) { buildSystemMatrix_impl(cavity); }
    /*! \brief Initialize CGSolver_ data member
     *  \param[in] max_it maximum number of iterations in the conjugate gradient solver
     *  \param[in] tol    tolerance in the conjugate gradient solver
     */
    void initializeCGSolver(int max_it, double tol) {
        if (cgBuilt_) PCMSOLVER_ERROR("Conjugate gradient solver already initialized", BOOST_CURRENT_FUNCTION);
        CGSolver_.setMaxIterations(max_it);
        CGSolver_.setTolerance(tol);
        CGSolver_.compute(TRdagger_);
        cgBuilt_ = true;
    }
    /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
     *  \param[in] dt the time step for the Euler integrator
     *  \param[in] MEP_current the vector containing the MEP at cavity points, at time (t + dt)
     *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time t
     *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time t
     */
    Eigen::VectorXd propagateASC(double dt, const Eigen::VectorXd & MEP_current, const Eigen::VectorXd & MEP_previous,
            const Eigen::VectorXd & ASC_previous) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
        if (!cgBuilt_) PCMSOLVER_ERROR("Conjugate gradient solver not initialized yet", BOOST_CURRENT_FUNCTION);
        return propagateASC_impl(dt, MEP_current, MEP_previous, ASC_previous);
    }
    /*! \brief Returns the ASC at initial time
     *  \param[in] MEP the vector containing the MEP at cavity points at initial time
     *  Uses dynamic PCM matrix to compute the initial value for the ASC
     */
    Eigen::VectorXd initialValueASC(const Eigen::VectorXd & MEP) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
        if (!cgBuilt_) PCMSOLVER_ERROR("Conjugate gradient solver not initialized yet", BOOST_CURRENT_FUNCTION);
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

    friend std::ostream & operator<<(std::ostream & os, TDPCMIterativeSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    /*! Time-dependent isotropic, uniform Debye permittivity profile */
    Debye permittivity_;
    /*! Whether the system matrix has been built */
    bool built_;
    /*! Whether the CG solver has been initialized */
    bool cgBuilt_;
    /*! PCM matrix */
    Eigen::MatrixXd TRdagger_;
    /*! Conjugate Gradient solver */
    Eigen::ConjugateGradient<Eigen::MatrixXd> CGSolver_;

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

#endif // TDPCMITERATIVESOLVER_HPP
