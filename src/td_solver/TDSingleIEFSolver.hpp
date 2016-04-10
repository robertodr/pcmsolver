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

#ifndef TDSINGLEIEFSOLVER_HPP
#define TDSINGLEIEFSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"
#include "Debye.hpp"
#include "TDPCMSolver.hpp"

/*! \file TDSingleIEFSolver.hpp
 *  \class TDSingleIEFSolver
 *  \brief Time-dependent solver for isotropic IEF with a single relaxation time
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class TDSingleIEFSolver : public TDPCMSolver
{
public:
    TDSingleIEFSolver() {}
    /*! \brief Construct solver from two Green's functions
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     *  \param[in] tau IEF relaxation time
     */
    TDSingleIEFSolver(double es, double ed, double t, double tau)
            : TDPCMSolver(es, ed, t), tauIEF_(tau) {}
    virtual ~TDSingleIEFSolver() {}
    friend std::ostream & operator<<(std::ostream & os, TDSingleIEFSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    /*! IEF relaxation time */
    double tauIEF_;
private:
    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity) __override {
        // The total size of the cavity
        int cavitySize = cavity.size();
        // Identiy matrix
        Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
        double e_d = permittivity_.epsilonDynamic;
        double e_0 = permittivity_.epsilonStatic;

        Vacuum<DerivativeTraits, IntegratorPolicy> gf_i;
        // Compute S and D on the whole cavity, regardless of symmetry
        Eigen::MatrixXd S = gf_i.singleLayer(cavity.elements());
        Eigen::MatrixXd D = gf_i.doubleLayer(cavity.elements());

        Eigen::MatrixXd A = cavity.elementArea().asDiagonal();

        // Form A_ (the dynamic matrix)
        double f_d = (e_d + 1.0) / (e_d - 1.0);
        A_ = 2 * M_PI * f_d * S - D * A * S;
        Eigen::FullPivLU<Eigen::MatrixXd> A__LU(A_);
        if (!(A__LU.isInvertible())) PCMSOLVER_ERROR("A_ matrix is not invertible!", BOOST_CURRENT_FUNCTION);
        A_ = -A__LU.inverse();
        A_ *= (2 * M_PI * Id - D * A);
        hermitivitize(A_);
        // Form B_ (the static matrix)
        double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
        B_ = 2 * M_PI * f_0 * S - D * A * S;
        Eigen::FullPivLU<Eigen::MatrixXd> B__LU(B_);
        if (!(B__LU.isInvertible())) PCMSOLVER_ERROR("B_ matrix is not invertible!", BOOST_CURRENT_FUNCTION);
        B_ = -B__LU.inverse();
        B_ *= (2 * M_PI * Id - D * A);
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
        double factor = dt / tauIEF_;
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
        os << "Solver Type: IEFPCM, isotropic" << std::endl;
        os << "IEF relaxation time = " << tauIEF_;
        return os;
    }
};

#endif // TDSINGLEIEFSOLVER_HPP
