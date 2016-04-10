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

#ifndef TDCPCMITERATIVESOLVER_HPP
#define TDCPCMITERATIVESOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"
#include "Debye.hpp"
#include "TDPCMIterativeSolver.hpp"


/*! \file TDCPCMIterativeSolver.hpp
 *  \class TDCPCMIterativeSolver
 *  \brief Time-dependent solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class TDCPCMIterativeSolver : public TDPCMIterativeSolver
{
public:
    TDCPCMIterativeSolver() {}
    /*! \brief Construct solver from two Green's functions
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     *  \param[in] corr factor to correct the conductor results
     */
    TDCPCMIterativeSolver(double es, double ed, double t, double corr)
            : TDPCMIterativeSolver(es, ed, t), correction_(corr) {}
    virtual ~TDCPCMIterativeSolver() {}
    friend std::ostream & operator<<(std::ostream & os, TDCPCMIterativeSolver & solver) {
        return solver.printSolver(os);
    }
private:
    /*! Correction for the conductor results */
    double correction_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity) __override {
        Vacuum<DerivativeTraits, IntegratorPolicy> gf_i;
        // Compute SI on the whole cavity, regardless of symmetry
        Eigen::MatrixXd SI = gf_i.singleLayer(cavity.elements());
        hermitivitize(SI);
        TRdagger_ = SI;

        built_ = true;
    }
    /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
     *  \param[in] dt the time step for the Euler integrator
     *  \param[in] MEP_current the vector containing the MEP at cavity points, at time (t + dt)
     *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time t
     *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time t
     */
    virtual Eigen::VectorXd propagateASC_impl(double dt, const Eigen::VectorXd & MEP_current, const Eigen::VectorXd & MEP_previous,
            const Eigen::VectorXd & ASC_previous) const __override {
        // Form RHS
        double f_d = permittivity_.dynamicOnsager(correction_);
        double f_0 = permittivity_.staticOnsager(correction_);
        double tau_CPCM = permittivity_.tau * (permittivity_.epsilonDynamic / permittivity_.epsilonStatic);
        double factor = dt / tau_CPCM;
        Eigen::VectorXd tilde_ASC = TRdagger_ * ASC_previous;
        Eigen::VectorXd RHS = - f_d * (MEP_current - MEP_previous) + tilde_ASC - factor * (tilde_ASC - f_0 * MEP_previous);

        Eigen::VectorXd ASC_current = CGSolver_.solve(RHS);

        return ASC_current;
    }
    /*! \brief Returns the ASC at initial time
     *  \param[in] MEP the vector containing the MEP at cavity points at initial time
     *  Uses dynamic PCM matrix to compute the initial value for the ASC
     */
    virtual Eigen::VectorXd initialValueASC_impl(const Eigen::VectorXd & MEP) const __override {
        double f_d = permittivity_.dynamicOnsager(correction_);

        return CGSolver_.solve(f_d * MEP);
    }
    virtual std::ostream & printSolver(std::ostream & os) {
        os << "Solver Type: iterative C-PCM" << std::endl;
        os << "Correction factor = " << correction_ << std::endl;
        os << "CG maximum iterations = " << CGSolver_.maxIterations() << std::endl;
        os << "CG tolerance = " << CGSolver_.tolerance();
        return os;
    }
};

#endif // TDCPCMITERATIVESOLVER_HPP
