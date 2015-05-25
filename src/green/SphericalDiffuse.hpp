/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
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

#ifndef SPHERICALDIFFUSE_HPP
#define SPHERICALDIFFUSE_HPP

#include <array>
#include <cmath>
#include <functional>
#include <iosfwd>
#include <string>
#include <tuple>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

// Has to be included here
#include "InterfacesImpl.hpp"
// Boost.Math includes
#include <boost/math/special_functions/legendre.hpp>

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "LoggerInterface.hpp"
#include "MathUtils.hpp"
#include "Timer.hpp"

/*! \file SphericalDiffuse.hpp
 *  \class SphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 *  \tparam ProfilePolicy functional form of the diffuse layer
 *
 *  This class is general, in the sense that no specific dielectric
 *  profile has been set in its definition.
 *  In principle any profile that can be described by:
 *  1. a left-side dielectric constant;
 *  2. a right-side dielectric constant;
 *  3. an interface layer width;
 *  4. an interface layer center
 *  can be used to define a new diffuse interface with spherical symmetry.
 */

template <typename IntegratorPolicy,
          typename ProfilePolicy>
class SphericalDiffuse : public GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy,
                                               SphericalDiffuse<IntegratorPolicy, ProfilePolicy> >
{
public:
    /*! Constructor for a one-layer interface
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     * \param[in] o center of the sphere
     */
    SphericalDiffuse(double e1, double e2, double w, double c, const Eigen::Vector3d & o)
        : GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, SphericalDiffuse<IntegratorPolicy, ProfilePolicy> >(), origin_(o)
    {
        initProfilePolicy(e1, e2, w, c);
        initSphericalDiffuse();
    }
    virtual ~SphericalDiffuse() {}
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator
     *  for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        double eps_r2 = 0.0;
        std::tie(eps_r2, std::ignore) = this->profile_(p2.norm());

        return (eps_r2 * this->derivativeProbe(direction, p1, p2));
    }

    /*! Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalS(const Element & e) const
    {
            return this->diagonal_.computeS(*this, e);
    }
    /*! Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalD(const Element & e) const
    {
            return this->diagonal_.computeD(*this, e);
    }

    friend std::ostream & operator<<(std::ostream & os, SphericalDiffuse & gf) {
        return gf.printObject(os);
    }
    /*! \brief Returns Coulomb singularity separation coefficient
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double coefficientCoulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        double r1  = (source - this->origin_).norm();
        double r2  = (probe  - this->origin_).norm();

        // Obtain coefficient for the separation of the Coulomb singularity
        return this->coefficient(r1, r2);
    }
    /*! \brief Returns singular part of the Green's function
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double Coulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        Eigen::Vector3d source_shifted = source - this->origin_;
        Eigen::Vector3d probe_shifted  = probe  - this->origin_;
        double r1  = source_shifted.norm();
        double r2  = probe_shifted.norm();
        double r12 = (source_shifted - probe_shifted).norm();

        // Obtain coefficient for the separation of the Coulomb singularity
        return (1.0 / (this->coefficient(r1, r2) * r12));
    }
    /*! \brief Returns non-singular part of the Green's function (image potential)
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double imagePotential(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        Eigen::Vector3d source_shifted = source - this->origin_;
        Eigen::Vector3d probe_shifted  = probe  - this->origin_;
        double r1  = source_shifted.norm();
        double r2  = probe_shifted.norm();
        double cos_gamma = source_shifted.dot(probe_shifted) / (r1 * r2);

        // Obtain coefficient for the separation of the Coulomb singularity
        double Cr12 = this->coefficient(r1, r2);

        double gr12 = 0.0;
        for (int L = 0; L <= maxLGreen_; ++L) {
            gr12 += this->functionSummation(L, r1, r2, cos_gamma, Cr12);
        }

        return gr12;
    }
    /*! Returns value of the directional derivative of the
     *  Coulomb singularity separation coefficient for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double coefficientCoulombDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        using namespace std::placeholders;
        return threePointStencil(std::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::coefficientCoulomb, this, _1, _2),
                p2, p1, direction, this->delta_);
    }
    /*! Returns value of the directional derivative of the
     *  singular part of the Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double CoulombDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        using namespace std::placeholders;
        return threePointStencil(std::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::Coulomb, this, _1, _2),
                p2, p1, direction, this->delta_);
    }
    /*! Returns value of the directional derivative of the
     *  non-singular part (image potential) of the Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double imagePotentialDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        using namespace std::placeholders;
        return threePointStencil(std::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::imagePotential, this, _1, _2),
                p2, p1, direction, this->delta_);
    }
    /*! Handle to the dielectric profile evaluation */
    void epsilon(double & v, double & d, double point) const { std::tie(v, d) = this->profile_(point); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    using RadialFunction = interfaces::RadialFunction;
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     *
     *  \note This takes care of the origin shift
     */
    virtual Numerical operator()(Numerical * sp, Numerical * pp) const
    {
        // Transfer raw arrays to Eigen vectors using the Map type
        Eigen::Map<Eigen::Matrix<double, 3, 1> > p1(sp), p2(pp);
        Eigen::Vector3d source = p1 - this->origin_;
        Eigen::Vector3d probe  = p2 - this->origin_;
        double r1  = source.norm();
        double r2  = probe.norm();
        double r12 = (source - probe).norm();
        double cos_gamma = source.dot(probe) / (r1 * r2);

        // Obtain coefficient for the separation of the Coulomb singularity
        double Cr12 = this->coefficient(r1, r2);

        double gr12 = 0.0;
        for (int L = 0; L <= maxLGreen_; ++L) {
            gr12 += this->functionSummation(L, r1, r2, cos_gamma, Cr12);
        }

        return (1.0 / (Cr12 * r12) + gr12);
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: spherical diffuse" << std::endl;
        os << this->profile_ << std::endl;
        os << "Sphere center = " << this->origin_.transpose() << std::endl;
        os << "Angular momentum (Green's function)    = " << this->maxLGreen_ << std::endl;
        os << "Angular momentum (Coulomb coefficient) = " << this->maxLC_;
        return os;
    }
    /*! Initializes a one-layer profile
     *  \param[in] e1 left-side dielectric constant
     *  \param[in] e2 right-side dielectric constant
     *  \param[in] w width of the interface layer
     *  \param[in] c center of the diffuse layer
     */
    void initProfilePolicy(double e1, double e2, double w, double c)
    { this->profile_ = ProfilePolicy(e1, e2, w, c); }
    /*! This calculates all the components needed to evaluate the Green's function */
    void initSphericalDiffuse() {
        using namespace std::placeholders;
        using namespace interfaces;

        LOG("SphericalDiffuse::initSphericalDiffuse");
        // Parameters for the numerical solution of the radial differential equation
        double r_0_         = 0.5;     /*! Lower bound of the integration interval */
        double r_infinity_  = this->profile_.center() + 200.0; /*! Upper bound of the integration interval */
        double observer_step_ = 1.0e-03; /*! Time step between observer calls */
        IntegratorParameters params_(r_0_, r_infinity_, observer_step_);
        ProfileEvaluator eval_ = std::bind(&TanhDiffuse::operator(), this->profile_, _1);

        LOG("Computing coefficient for the separation of the Coulomb singularity");
        LOG("Computing first radial solution L = " + std::to_string(maxLC_));
        timerON("computeZeta for coefficient");
        computeZeta(maxLC_, zetaC_, eval_, params_);
        timerOFF("computeZeta for coefficient");
        LOG("DONE: Computing first radial solution L = " + std::to_string(maxLC_));

        LOG("Computing second radial solution L = " + std::to_string(maxLC_));
        timerON("computeOmega for coefficient");
        computeOmega(maxLC_, omegaC_, eval_, params_);
        timerOFF("computeOmega for coefficient");
        LOG("Computing second radial solution L = " + std::to_string(maxLC_));
        LOG("DONE: Computing coefficient for the separation of the Coulomb singularity");

        LOG("Computing radial solutions for Green's function");
        timerON("SphericalDiffuse: Looping over angular momentum");
        for (int L = 0; L <= maxLGreen_; ++L) {
            // First radial solution
            LOG("Computing first radial solution L = " + std::to_string(L));
            timerON("computeZeta L = " + std::to_string(L));
            // Create an empty RadialFunction
            RadialFunction tmp_zeta_;
            computeZeta(L, tmp_zeta_, eval_, params_);
            zeta_.push_back(tmp_zeta_);
            timerOFF("computeZeta L = " + std::to_string(L));
            LOG("DONE: Computing first radial solution L = " + std::to_string(L));

            // Second radial solution
            LOG("Computing second radial solution L = " + std::to_string(L));
            timerON("computeOmega L = " + std::to_string(L));
            // Create an empty RadialFunction
            RadialFunction tmp_omega_;
            computeOmega(L, tmp_omega_, eval_, params_);
            omega_.push_back(tmp_omega_);
            timerOFF("computeOmega L = " + std::to_string(L));
            LOG("DONE: Computing second radial solution L = " + std::to_string(L));
        }
        timerOFF("SphericalDiffuse: Looping over angular momentum");
        LOG("DONE: Computing radial solutions for Green's function");
    }

    /*! Center of the dielectric sphere */
    Eigen::Vector3d origin_;

    /**@{ Parameters and functions for the calculation of the Green's function, including Coulomb singularity */
    /*! Maximum angular momentum in the final summation over Legendre polynomials to obtain G */
    int maxLGreen_ = 30;
    /*! \brief First independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_ and has r^l behavior
     */
    std::vector<RadialFunction> zeta_;
    /*! \brief Second independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_  and has r^(-l-1) behavior
     */
    std::vector<RadialFunction> omega_;
    /*! \brief Returns L-th component of the radial part of the Green's function
     *  \param[in] L  angular momentum
     *  \param[in] r1 first point
     *  \param[in] r2 second point
     *  \param[in] cos_gamma cosine of the angle in between first and second point
     *  \param[in] Cr12 Coulomb singularity separation coefficient
     */
    double functionSummation(int L, double r1, double r2, double cos_gamma, double Cr12) const {
        double gr12 = 0.0;
        // Evaluate Legendre polynomial of order L
        // First of all clean-up cos_gamma, Legendre polynomials
        // are only defined for -1 <= x <= 1
        if (numericalZero(cos_gamma - 1)) cos_gamma = 1.0;
        if (numericalZero(cos_gamma + 1)) cos_gamma = -1.0;
        double pl_x = boost::math::legendre_p(L, cos_gamma);

        /* Value of zeta_[L] at point with index 1 */
        double zeta1  = linearInterpolation(r1, zeta_[L][0], zeta_[L][1]);
        /* Value of zeta_[L} at point with index 2 */
        double zeta2  = linearInterpolation(r2, zeta_[L][0], zeta_[L][1]);
        /* Value of omega_[L] at point with index 1 */
        double omega1 = linearInterpolation(r1, omega_[L][0], omega_[L][1]);
        /* Value of omega_[L} at point with index 2 */
        double omega2 = linearInterpolation(r2, omega_[L][0], omega_[L][1]);

        /* Components for the evaluation of the Wronskian */
        /* Value of derivative of zeta_[L] at point with index 2 */
        double d_zeta2  = linearInterpolation(r2, zeta_[L][0], zeta_[L][2]);
        /* Value of derivative of omega_[L] at point with index 2 */
        double d_omega2 = linearInterpolation(r2, omega_[L][0], omega_[L][2]);

        double eps_r2 = 0.0;
        std::tie(eps_r2, std::ignore) = this->profile_(r2);

        double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

        if (r1 < r2) {
            gr12 = std::exp(zeta1 - zeta2) * (2*L +1) / denominator;
            gr12 = (gr12 - std::pow(r1/r2, L) / (r2 * Cr12) ) * pl_x ;
        } else {
            gr12 = std::exp(omega1 - omega2) * (2*L +1) / denominator;
            gr12 = (gr12 - std::pow(r2/r1, L) / (r1 * Cr12) ) * pl_x ;
        }

        return gr12;
    }
    /**@}*/

    /**@{ Parameters and functions for the calculation of the Coulomb singularity separation coefficient */
    /*! Maximum angular momentum to obtain C(r, r'), needed to separate the Coulomb singularity */
    int maxLC_     = 2 * maxLGreen_;
    /*! \brief First independent radial solution, used to build coefficient.
     *  \note This is needed to separate the Coulomb singularity and has r^l behavior
     */
    RadialFunction zetaC_;
    /*! \brief Second independent radial solution, used to build coefficient.
     *  \note This is needed to separate the Coulomb singularity and has r^(-l-1) behavior
     */
    RadialFunction omegaC_;
    /*! \brief Returns coefficient for the separation of the Coulomb singularity
     *  \param[in] r1 first point
     *  \param[in] r2 second point
     */
    double coefficient(double r1, double r2) const {
        /* Value of zetaC_ at point with index 1 */
        double zeta1  = linearInterpolation(r1, zetaC_[0], zetaC_[1]);
        /* Value of zetaC_ at point with index 2 */
        double zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[1]);
        /* Value of omegaC_ at point with index 1 */
        double omega1 = linearInterpolation(r1, omegaC_[0], omegaC_[1]);
        /* Value of omegaC_ at point with index 2 */
        double omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[1]);

        /* Components for the evaluation of the Wronskian */
        /* Value of derivative of zetaC_ at point with index 2 */
        double d_zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[2]);
        /* Value of derivative of omegaC_ at point with index 2 */
        double d_omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[2]);

        double tmp = 0.0, coeff = 0.0;
        double eps_r2 = 0.0;
        std::tie(eps_r2, std::ignore) = this->profile_(r2);

        double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

        if (r1 < r2) {
            tmp = std::exp(zeta1 - zeta2) * (2*maxLC_ +1) / denominator;
            coeff = std::pow(r1/r2, maxLC_) / (tmp * r2);
        } else {
            tmp = std::exp(omega1 - omega2) * (2*maxLC_ +1) / denominator;
            coeff = std::pow(r2/r1, maxLC_) / (tmp * r1);
        }

        return coeff;
    }
    /**@}*/
};

#endif // SPHERICALDIFFUSE_HPP
