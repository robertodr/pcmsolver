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

#ifndef TDIEFSOLVER_HPP
#define TDIEFSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "Cavity.hpp"
#include "Debye.hpp"
#include "Element.hpp"
#include "MathUtils.hpp"
#include "Vacuum.hpp"
#include "TDSolverHelperFunctions.hpp"
#include "TDPCMSolver.hpp"

/*! \file TDIEFSolver.hpp
 *  \class TDIEFSolver
 *  \brief Time-dependent solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits,
          typename IntegratorPolicy>
class TDIEFSolver : public TDPCMSolver
{
public:
    TDIEFSolver() {}
    /*! \brief Construct solver from two Green's functions
     *  \param[in] es static permittivity
     *  \param[in] ed dynamic permittivity
     *  \param[in] t  relaxation time
     *  \param[in] cholesky whether Cholesky or Lowdin decomposition is used
     */
    TDIEFSolver(double es, double ed, double t, bool cholesky)
            : TDPCMSolver(es, ed, t), useCholesky_(cholesky) {}
    virtual ~TDIEFSolver() {}
    /*! Return i-th element of the static K matrix
     *  \param[in] i index
     */
    double K_0(int i) const { return K_0_(i); }
    /*! Return i-th element of the dynamic K matrix
     *  \param[in] i index
     */
    double K_d(int i) const { return K_d_(i); }
    /*! Return i-th element of the Lambda matrix
     *  \param[in] i index
     */
    double Lambda(int i) const { return Lambda_(i); }
    /*! Return i-th element of the relaxation times matrix
     *  \param[in] i index
     *  \note Relaxation times are expressed in atomic units
     */
    double tau(int i) const { return tau_(i); }
    /*! Return static K matrix */
    const Eigen::VectorXd & K_0() const { return K_0_; }
    /*! Return dynamic K matrix */
    const Eigen::VectorXd & K_d() const { return K_d_; }
    /*! Return Lambda matrix */
    const Eigen::VectorXd & Lambda() const { return Lambda_; }
    /*! Return relaxation times matrix
     *  \note Relaxation times are expressed in atomic units
     */
    const Eigen::VectorXd & tau() const { return tau_; }
    friend std::ostream & operator<<(std::ostream & os, TDIEFSolver & solver) {
        return solver.printSolver(os);
    }
private:
    /*! Whether Cholesky decomposition is used to solve the generalized eigenvalue problem */
    bool useCholesky_;
    /*! Diagonalized PCM matrix, static */
    Eigen::VectorXd K_0_;
    /*! Diagonalized PCM matrix, dynamic */
    Eigen::VectorXd K_d_;
    /*! Diagonalized S^-1/2DAS^1/2 matrix */
    Eigen::VectorXd Lambda_;
    /*! Diagonal relaxation times matrix */
    Eigen::VectorXd tau_;
    /*! Relaxation times (used in the propagator) */
    Eigen::MatrixXd C_;
    /*! Dynamic isotropic IEF matrix (used for the initial value of the ASC) */
    Eigen::MatrixXd PCMMatrix_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity) __override {
        if (useCholesky_) {
            systemMatrix_Cholesky(cavity);
        } else {
            systemMatrix_Lowdin(cavity);
        }
    }
    void systemMatrix_Lowdin(const Cavity & cavity) {
        using namespace td_solver;
        // The total size of the cavity
        int cavitySize = cavity.size();
        // Diagonal areas matrix
        Eigen::MatrixXd A = cavity.elementArea().asDiagonal();
        // Identiy matrix
        Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
        double e_d = permittivity_.epsilonDynamic;
        double e_0 = permittivity_.epsilonStatic;

        Vacuum<DerivativeTraits, IntegratorPolicy> gf_i;
        // Compute S on the whole cavity, regardless of symmetry
        Eigen::MatrixXd S = gf_i.singleLayer(cavity.elements());
        // Compute D on whole cavity, regardless of symmetry
        Eigen::MatrixXd D = gf_i.doubleLayer(cavity.elements());
        // Make sure S is symmetric (might be unnecessary)
        hermitivitize(S);
        // Compute S^1/2 and S^-1/2
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> S_eigen(S);
        Eigen::MatrixXd S_sqrt = S_eigen.operatorSqrt();
        Eigen::MatrixXd S_invsqrt = S_eigen.operatorInverseSqrt();
        // Compute S^-1/2DAS^1/2 and symmetrize it
        D = S_invsqrt * D * A * S_sqrt;
        hermitivitize(D);
        // Diagonalize S^-1/2DAS^1/2
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> D_eigen(D);
        // Get Lambda and the transformation matrix
        Lambda_ = D_eigen.eigenvalues();
        Eigen::MatrixXd T = D_eigen.eigenvectors();
        // Form tau^-1 (relaxation times matrix)
        tau_ = td_solver::tau(Lambda_, e_d, e_0, permittivity_.tau);
        Eigen::VectorXd tau_inv = tauInverse(Lambda_, e_d, e_0, permittivity_.tau);
        // Form A_ (the dynamic matrix)
        double f_d = (e_d + 1.0) / (e_d - 1.0);
        K_d_ = K(Lambda_, f_d);
        A_ = - S_invsqrt * T * K_d_.asDiagonal() * T.adjoint().eval() * S_invsqrt;
        // Form B_ (the stati matrix)
        double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
        K_0_ = K(Lambda_, f_0);
        B_ = - S_invsqrt * T * tau_inv.asDiagonal() * K_0_.asDiagonal() * T.adjoint().eval() * S_invsqrt;
        // Form C_ (the relaxation times matrix)
        C_ = S_invsqrt * T * tau_inv.asDiagonal() * T.adjoint().eval() * S_sqrt;
        // Form PCMMatrix_ (needed to form the initial value of the ASC)
        PCMMatrix_ = - S_invsqrt * T * K_d_.asDiagonal() * T.adjoint().eval() * S_invsqrt;

        built_ = true;
    }
    void systemMatrix_Cholesky(const Cavity & /* cavity */) {
        /*
        using namespace td_solver;
        // The total size of the cavity
        int cavitySize = cavity.size();
        // Diagonal areas matrix
        Eigen::MatrixXd A = cavity.elementArea().asDiagonal();
        // Identiy matrix
        Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
        double e_d = permittivity_.epsilonDynamic;
        double e_0 = permittivity_.epsilonStatic;

        Vacuum<DerivativeTraits, IntegratorPolicy> gf_i;
        // Compute S on the whole cavity, regardless of symmetry
        Eigen::MatrixXd S = gf_i.singleLayer(cavity.elements());
        // Compute D on whole cavity, regardless of symmetry
        Eigen::MatrixXd D = gf_i.doubleLayer(cavity.elements());
        // Make sure S is symmetric (might be unnecessary)
        hermitivitize(S);
        // Compute Cholesky decomposition of S = LVL^t
        Eigen::LDLT<Eigen::MatrixXd> S_lvlt(S);
        Eigen::MatrixXd L = S_lvlt.matrixL();
        Eigen::MatrixXd L_inv = L.inverse();
        Eigen::MatrixXd Lt_inv = L_inv.transpose();
        Eigen::MatrixXd V = S_lvlt.vectorD().asDiagonal();
        Eigen::MatrixXd V_inv = V.inverse();
        // Compute L^-1DAL and symmetrize it
        D = L.inverse() * D * A * L;
        //hermitivitize(D);
        // Diagonalize L^-1DAL
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> D_eigen(D);
        // Get Lambda and the transformation matrix
        Eigen::MatrixXd Lambda = D_eigen.eigenvalues().asDiagonal();
        std::cout << D_eigen.eigenvalues() << std::endl;
        Eigen::MatrixXd T = D_eigen.eigenvectors();
        // Form tau^-1 (relaxation times matrix)
        Eigen::MatrixXd tau_inv = tauInverse(Lambda, e_d, e_0, permittivity_.tau);
        // Form A_ (the dynamic matrix)
        double f_d = (e_d + 1.0) / (e_d - 1.0);
        A_ = - Lt_inv * T * K(Lambda, f_d) * T.adjoint().eval() * V_inv * L_inv;
        // Form B_ (the static matrix)
        double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
        B_ = - Lt_inv * T * tau_inv * K(Lambda, f_0) * T.adjoint().eval() * V_inv * L_inv;
        // Form C_ (the relaxation times matrix)
        C_ = Lt_inv * T * tau_inv * T.adjoint().eval() * L * V;
        // Form PCMMatrix_ (needed to form the initial value of the ASC)
        PCMMatrix_ = - Lt_inv * T * K(Lambda, f_d) * T.adjoint().eval() * V_inv * L_inv;

        built_ = true;
        */
    }
    /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
     *  \param[in] dt the time step for the Euler integrator
     *  \param[in] MEP_current the vector containing the MEP at cavity points, at time (t + dt)
     *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time t
     *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time t
     */
    virtual Eigen::VectorXd propagateASC_impl(double dt, const Eigen::VectorXd & MEP_current, const Eigen::VectorXd & MEP_previous,
            const Eigen::VectorXd & ASC_previous) const __override {
        return (A_ * (MEP_current - MEP_previous) + dt * B_ * MEP_previous - dt * C_ * ASC_previous + ASC_previous);
    }
    /*! \brief Returns the ASC at initial time
     *  \param[in] MEP the vector containing the MEP at cavity points at initial time
     *  Uses dynamic PCM matrix to compute the initial value for the ASC
     */
    virtual Eigen::VectorXd initialValueASC_impl(const Eigen::VectorXd & MEP) const __override {
        return (PCMMatrix_ * MEP);
    }
    virtual std::ostream & printSolver(std::ostream & os) __override {
        os << "Solver Type: IEFPCM, isotropic" << std::endl;
        os << permittivity_;
        return os;
    }
};

#endif // TDIEFSOLVER_HPP
