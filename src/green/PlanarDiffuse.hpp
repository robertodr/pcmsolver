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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#ifndef PLANARDIFFUSE_HPP
#define PLANARDIFFUSE_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

// Has to be included here
#include "InterfacesImpl.hpp"
// Boost.Math includes
#include <boost/math/special_functions/bessel.hpp>

#include "DerivativeUtils.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "dielectric_profile/ProfileForward.hpp"
#include "GreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "utils/QuadratureRules.hpp"

/*! \file PlanarDiffuse.hpp
 *  \class PlanarDiffuse
 *  \brief Green's function for a diffuse interface with planar symmetry
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2016
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 *  \tparam ProfilePolicy functional form of the diffuse layer
 *
 *  This class is general, in the sense that no specific dielectric profile has been set in its definition.
 *  In principle any profile that can be described by:
 *  1. a left-side dielectric constant;
 *  2. a right-side dielectric constant;
 *  3. an interface layer width;
 *  4. an interface layer center
 *  can be used to define a new diffuse interface with cylindrical symmetry,
 *  i.e. an infinite plane.
 *  The origin of the dielectric plane can be changed by means of the constructor.
 *  The solution of the differential equation defining the Green's function is **always**
 *  performed assuming that the dielectric plane is centered in the origin of the coordinate
 *  system. Whenever the public methods are invoked to "sample" the Green's function
 *  at a pair of points, a translation of the sampling points is performed first.
 */

template <typename IntegratorPolicy = CollocationIntegrator,
          typename ProfilePolicy = OneLayerTanh>
class PlanarDiffuse __final : public GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy,
                                               PlanarDiffuse<IntegratorPolicy, ProfilePolicy> >
{
public:
    /*! Constructor for a one-layer interface
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] zInt position of the diffuse layer
     */
    PlanarDiffuse(double e1, double e2, double w, double zInt)
        : GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, PlanarDiffuse<IntegratorPolicy, ProfilePolicy> >(),
		origin_(Eigen::Vector3d::UnitZ() * zInt)
    {
        initProfilePolicy(e1, e2, w, zInt);
        initPlanarDiffuse();
    }
    /*! Constructor for a one-layer interface
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] zInt position of the diffuse layer
     * \param[in] f scaling factor for diagnal elements
     */
    PlanarDiffuse(double e1, double e2, double w, double zInt, double f)
        : GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, PlanarDiffuse<IntegratorPolicy, ProfilePolicy> >(f),
		origin_(Eigen::Vector3d::UnitZ() * zInt)
    {
        initProfilePolicy(e1, e2, w, zInt);
        initPlanarDiffuse();
    }
    virtual ~PlanarDiffuse() {}

    /*! Calculates the matrix representation of the S operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const __override
    {
        return this->integrator_.singleLayer(*this, e);
    }
    /*! Calculates the matrix representation of the D operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const __override
    {
        return this->integrator_.doubleLayer(*this, e);
    }

    friend std::ostream & operator<<(std::ostream & os, PlanarDiffuse & gf) {
        return gf.printObject(os);
    }
    /*! \brief Returns Coulomb singularity separation coefficient
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double coefficientCoulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        // Obtain coefficient for the separation of the Coulomb singularity
        return this->coefficient_impl(source, probe);
    }
    /*! \brief Returns singular part of the Green's function
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double Coulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        double r12 = (source - probe).norm();

        // Obtain coefficient for the separation of the Coulomb singularity
        return (1.0 / (this->coefficient_impl(source, probe) * r12));
    }
    /*! \brief Returns non-singular part of the Green's function (image potential)
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double imagePotential(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        // Obtain coefficient for the separation of the Coulomb singularity
        double Cr12 = this->coefficient_impl(source, probe);
		double rho = ((source-probe).head(2)).norm();
		double greenImage = 0.0;
		
        for (int kindex = 0; kindex < 64; kindex++) {
			double gauss_point = rule_.gaussAbscissa(kindex);
			double kstep = -std::log((gauss_point + 1.0)/2.0);
			double bess_0_x = boost::math::cyl_bessel_j(0, kstep * rho);
            greenImage += kstep * bess_0_x * this->imagePotentialComponent_impl(kindex, source, probe, Cr12);
        }

        return greenImage;

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
        return threePointStencil(pcm::bind(&PlanarDiffuse<IntegratorPolicy, ProfilePolicy>::coefficientCoulomb, this, pcm::_1, pcm::_2),
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
        return threePointStencil(pcm::bind(&PlanarDiffuse<IntegratorPolicy, ProfilePolicy>::Coulomb, this, pcm::_1, pcm::_2),
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
        return threePointStencil(pcm::bind(&PlanarDiffuse<IntegratorPolicy, ProfilePolicy>::imagePotential, this, pcm::_1, pcm::_2),
                p2, p1, direction, this->delta_);
    }
    /*! Handle to the dielectric profile evaluation */
    pcm::tuple<double, double> epsilon(const Eigen::Vector3d & point) const {
        return this->profile_((point + this->origin_).norm());
    }
    void toFile(const std::string & prefix = "") {
		std::string tmp;
		prefix.empty() ? tmp = prefix : tmp = prefix + "-";
		for (int kindex = 0; kindex < 64; kindex++) {
			writeToFile(U1_[kindex], tmp + "U1_" + pcm::to_string(kindex) + ".dat");
			writeToFile(U2_[kindex], tmp + "U2_" + pcm::to_string(kindex) + ".dat");
		}
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     *
     *  \note This takes care of the origin shift
     */
    virtual Numerical operator()(Numerical * sp, Numerical * pp) const __override
    {
        // Transfer raw arrays to Eigen vectors using the Map type
        Eigen::Map<Eigen::Matrix<double, 3, 1> > source(sp), probe(pp);
        // Obtain coefficient for the separation of the Coulomb singularity
		double r12 = (source - probe).norm();
        double Cr12 = this->coefficient_impl(source, probe);
        double gr12 = this->imagePotential(source, probe);
		std::cout << "r12 Cr12 gr12 " << r12<< " " << Cr12 << " " << gr12 << std::endl;
        return (1.0 / (Cr12 * r12) + gr12);
    }

	gauss_bilateral64 rule_;
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __override
    {
        double eps_r2 = 0.0;
        // Shift p2 by origin_
        pcm::tie(eps_r2, pcm::ignore) = this->epsilon(p2);

        return (eps_r2 * this->derivativeProbe(direction, p1, p2));
    }
    virtual KernelS exportKernelS_impl() const __override {
      return pcm::bind(&PlanarDiffuse<IntegratorPolicy, ProfilePolicy>::kernelS, *this, pcm::_1, pcm::_2);
    }
    virtual KernelD exportKernelD_impl() const __override {
      return pcm::bind(&PlanarDiffuse<IntegratorPolicy, ProfilePolicy>::kernelD, *this, pcm::_1, pcm::_2, pcm::_3);
    }
    virtual std::ostream & printObject(std::ostream & os) __override
    {
        Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "(", ")");
        os << "Green's function type: planar diffuse" << std::endl;
        os << this->profile_ << std::endl;
        os << "Plane position        = " << this->origin_.transpose().format(CleanFmt) << std::endl;
        return os;
    }
    /*! Initializes a one-layer profile
     *  \param[in] e1 left-side dielectric constant
     *  \param[in] e2 right-side dielectric constant
     *  \param[in] w width of the interface layer
     *  \param[in] c center of the diffuse layer
     */
    void initProfilePolicy(double e1, double e2, double w, double c) {
        this->profile_ = ProfilePolicy(e1, e2, w, c);
    }
    /*! This calculates all the components needed to evaluate the Green's function */
    void initPlanarDiffuse() {
        using namespace interfaces;

        LOG("PlanarDiffuse::initPlanarDiffuse");
        // Parameters for the numerical solution of the radial differential equation
        double eps_abs_     = 1.0e-10; /*! Absolute tolerance level */
        double eps_rel_     = 1.0e-06; /*! Relative tolerance level */
        double factor_x_    = 0.0;     /*! Weight of the state      */
        double factor_dxdt_ = 0.0;     /*! Weight of the state derivative */
        double zmin_        = - 30.0; /*! Lower bound of the integration interval */
        double zmax_        =   30.0; /*! Upper bound of the integration interval */
        double observer_step_ = 1.0e-03; /*! Time step between observer calls */
        IntegratorParameters params_(eps_abs_, eps_rel_, factor_x_, factor_dxdt_, zmin_, zmax_, observer_step_);
        ProfileEvaluator eval_ = pcm::bind(&ProfilePolicy::operator(), this->profile_, pcm::_1);

        LOG("Computing radial solutions for Green's function");
        TIMER_ON("PlanarDiffuse: Looping over angular momentum");
        U1_.reserve(64);
        U2_.reserve(64);
        for (int kindex = 0; kindex < 64; kindex++) {
			std::cout << "kindex " << kindex << std::endl; 
			double gauss_point = rule_.gaussAbscissa(kindex);
			double kstep = -std::log((gauss_point + 1.0)/2.0);
            // First radial solution
            LOG("Computing first normal solution k = " + pcm::to_string(kstep));
			std::cout << "Computing first normal solution k = " << pcm::to_string(kstep) << std::endl;
            TIMER_ON("computeNormal 1 k = " + pcm::to_string(kstep));
            // Create an empty NormalFunction
            NormalFunction<StateType, NormalDifferential> tmp1_(kstep, zmin_, zmax_, eval_, params_);
			writeToFile(tmp1_, "U1_" + pcm::to_string(kindex) + ".dat");
            U1_.push_back(tmp1_);
            TIMER_OFF("computeNormal 1 k = " + pcm::to_string(kstep));
            LOG("DONE: Computing first normal solution k = " + pcm::to_string(kstep));

            // Second radial solution
            LOG("Computing second normal solution k = " + pcm::to_string(kstep));
			std::cout << "Computing second normal solution k = " << pcm::to_string(kstep) << std::endl;
            TIMER_ON("computeNormal 2 k = " + pcm::to_string(kstep));
            // Create an empty NormalFunction
            NormalFunction<StateType, NormalDifferential> tmp2_(kstep, zmax_, zmin_, eval_, params_);
			writeToFile(tmp2_, "U2_" + pcm::to_string(kindex) + ".dat");
            U2_.push_back(tmp2_);
            TIMER_OFF("computeNormal 2 k = " + pcm::to_string(kstep));
            LOG("DONE: Computing second normal solution k = " + pcm::to_string(kstep));
        }
        TIMER_OFF("PlanarDiffuse: Looping over angular momentum");
        LOG("DONE: Computing radial solutions for Green's function");
    }

    Eigen::Vector3d origin_;

    /**@{ Parameters and functions for the calculation of the Green's function, including Coulomb singularity */
    /*! \brief First independent normal solution, used to build Green's function.
     */
    std::vector<NormalFunction<interfaces::StateType, interfaces::NormalDifferential> > U1_;
    /*! \brief Second independent normal solution, used to build Green's function.
     */
    std::vector<NormalFunction<interfaces::StateType, interfaces::NormalDifferential> > U2_;
    /*! \brief Returns L-th component of the radial part of the Green's function
     *  \param[in] kindex index of the integration gridpoint
     *  \param[in] sp source point
     *  \param[in] pp probe point
     *  \param[in] Cr12 Coulomb singularity separation coefficient
     *  \note This function shifts the given source and probe points by the location of the
     *  dielectric plane.
     */
    double imagePotentialComponent_impl(int kindex, const Eigen::Vector3d & sp, const Eigen::Vector3d & pp, double Cr12) const {
		double gauss_point = rule_.gaussAbscissa(kindex);
		double kstep = -std::log((gauss_point + 1.0)/2.0);
		return G_k(kindex, sp, pp) - std::exp(-kstep*std::abs(sp(2)-pp(2))) / (Cr12 * kstep);
    }
    /**@}*/

    /**@{ Parameters and functions for the calculation of the Coulomb singularity separation coefficient */
    /*! \brief Returns coefficient for the separation of the Coulomb singularity
     *  \param[in] sp first point
     *  \param[in] pp second point
     *  \note This function shifts the given source and probe points by the location of the
     *  dielectric plane.
     */
    double coefficient_impl(const Eigen::Vector3d & sp, const Eigen::Vector3d & pp) const {
		double gauss_point = rule_.gaussAbscissa(63);
		double kstep = -std::log((gauss_point + 1.0)/2.0);
		double tmp = kstep * std::exp(kstep*std::abs(sp(2)-pp(2))) * G_k(63, sp, pp);
		std::cout << "gauss_point " << gauss_point << std::endl;
		std::cout << "kstep " << kstep << std::endl;
		std::cout << "tmp " << tmp << std::endl;
		for (int i = 0; i < 64; i++) {
			double GK =  G_k(i, sp, pp);
			std::cout << "GK " << i  << " " << GK << std::endl;
		}
		return 1.0/tmp;
    }
    /**@}*/

    double G_k(int kindex, const Eigen::Vector3d & sp, const Eigen::Vector3d & pp) const {
        Eigen::Vector3d sp_shift = sp + this->origin_;
        Eigen::Vector3d pp_shift = pp + this->origin_;

        /* Sample U1_ */
        double U1s = 0.0, U1p = 0.0, dU1p = 0.0;
        /* Value of zetaC_ at point with index 1 */
        pcm::tie(U1s, pcm::ignore) = U1_[kindex](sp_shift(2));
        /* Value of zetaC_ and its first derivative at point with index 2 */
        pcm::tie(U1p, dU1p) = U1_[kindex](pp_shift(2));

        /* Sample U2_ */
        double U2s = 0.0, U2p = 0.0, dU2p = 0.0;
        /* Value of zetaC_ at point with index 1 */
        pcm::tie(U2s, pcm::ignore) = U2_[kindex](sp_shift(2));
        /* Value of zetaC_ and its first derivative at point with index 2 */
        pcm::tie(U2p, dU2p) = U2_[kindex](pp_shift(2));

        double eps_r2 = 0.0;
        pcm::tie(eps_r2, pcm::ignore) = this->profile_(pp_shift(2));

        /* Evaluation of the Wronskian and the denominator */
        double denominator = eps_r2 * (dU2p - dU1p);

		double G_k = 0.0;

		if (pp_shift(2) < sp_shift(2)) {
			G_k = 2.0 * std::exp(U1s - U1p) / denominator;
        } else {
			G_k = 2.0 * std::exp(U2s - U2p) / denominator;
        }

        return G_k;
		
	}
};

#endif // PLANARDIFFUSE_HPP


