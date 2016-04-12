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

#ifndef TDCPCMSOLVER_HPP
#define TDCPCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "bi_operators/IntegratorTypes.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"
#include "Debye.hpp"
#include "TDPCMSolver.hpp"

/*! \file TDCPCMSolver.hpp
 *  \class TDCPCMSolver
 *  \brief Time-dependent solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class TDCPCMSolver : public TDPCMSolver
{
public:
    TDCPCMSolver() {}
    /*! \brief Construct solver from two Green's functions
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     *  \param[in] corr factor to correct the conductor results
     */
    TDCPCMSolver(double es, double ed, double t, double corr)
            : TDPCMSolver(es, ed, t), correction_(corr) {}
    virtual ~TDCPCMSolver() {}
    friend std::ostream & operator<<(std::ostream & os, TDCPCMSolver & solver) {
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

        // Invert SI  using LU decomposition with full pivoting
        // This is a rank-revealing LU decomposition, this allows us
        // to test if SI is invertible before attempting to invert it.
        Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
        if (!(SI_LU.isInvertible())) PCMSOLVER_ERROR("SI matrix is not invertible!", BOOST_CURRENT_FUNCTION);
        double f_d = permittivity_.dynamicOnsager(correction_);
        A_ = f_d * SI_LU.inverse();
        hermitivitize(A_);
        double f_0 = permittivity_.staticOnsager(correction_);
        B_ = f_0 * SI_LU.inverse();
        hermitivitize(B_);

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
        double tau_CPCM = permittivity_.tau * (permittivity_.epsilonDynamic / permittivity_.epsilonStatic);
        double factor = dt / tau_CPCM;
        return (A_ * (MEP_current - MEP_previous) + factor * B_ * MEP_previous - factor * ASC_previous + ASC_previous);
    }
    /*! \brief Returns the ASC at initial time
     *  \param[in] MEP the vector containing the MEP at cavity points at initial time
     *  Uses dynamic PCM matrix to compute the initial value for the ASC
     */
    virtual Eigen::VectorXd initialValueASC_impl(const Eigen::VectorXd & MEP) const __override {
        return (A_ * MEP);
    }
    virtual std::ostream & printSolver(std::ostream & os) __override {
        os << "Solver Type: C-PCM" << std::endl;
        os << "Correction factor = " << correction_;
        return os;
    }
};

#endif // TDCPCMSOLVER_HPP
