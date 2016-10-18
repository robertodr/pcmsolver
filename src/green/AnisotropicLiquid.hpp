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

#ifndef ANISOTROPICLIQUID_HPP
#define ANISOTROPICLIQUID_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "dielectric_profile/Anisotropic.hpp"
#include "GreensFunction.hpp"

/*! \file AnisotropicLiquid.hpp
 *  \class AnisotropicLiquid
 *  \brief Green's functions for anisotropic liquid, described by a tensorial
 * permittivity
 *  \author Roberto Di Remigio
 *  \date 2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */
// TODO Check the implementation of the coefficient_impl function!

template <typename DerivativeTraits = AD_directional>
class AnisotropicLiquid __final
    : public GreensFunction<DerivativeTraits, Anisotropic> {
public:
  /*! \param[in] eigen_eps eigenvalues of the permittivity tensors
   *  \param[in] euler_ang Euler angles in degrees
   */
  AnisotropicLiquid(const Eigen::Vector3d & eigen_eps,
                    const Eigen::Vector3d & euler_ang)
      : GreensFunction<DerivativeTraits, Anisotropic>() {
    this->profile_ = Anisotropic(eigen_eps, euler_ang);
  }
  virtual ~AnisotropicLiquid() {}

  friend std::ostream & operator<<(std::ostream & os, AnisotropicLiquid & gf) {
    return gf.printObject(os);
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See
                                     http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
                                     */
      private : virtual DerivativeTraits
                operator()(DerivativeTraits * source,
                           DerivativeTraits * probe) const __override {
    // The distance has to be calculated using epsilonInv_ as metric:
    DerivativeTraits scratch = 0.0;
    Eigen::Matrix3d epsilonInv_ = this->profile_.epsilonInv();
    double detEps_ = this->profile_.detEps();
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        scratch +=
            (source[i] - probe[i]) * epsilonInv_(i, j) * (source[j] - probe[j]);
      }
    }
    DerivativeTraits distance = sqrt(scratch);

    return (1.0 / (sqrt(detEps_) * distance));
  }
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const __override {
    // Since the permittivity is a tensorial quantity,
    // the full gradient is needed to get the kernel of D and D^\dagger
    Eigen::Vector3d scratch =
        this->profile_.epsilon() * (this->gradientProbe(p1, p2));
    return scratch.dot(direction);
  }

  virtual KernelS exportKernelS_impl() const __override {
    return pcm::bind(&AnisotropicLiquid<DerivativeTraits>::kernelS, *this, pcm::_1,
                     pcm::_2);
  }
  virtual KernelD exportKernelD_impl() const __override {
    return pcm::bind(&AnisotropicLiquid<DerivativeTraits>::kernelD, *this, pcm::_1,
                     pcm::_2, pcm::_3);
  }

  virtual double coefficient_impl(const Eigen::Vector3d & UNUSED(source),
                                  const Eigen::Vector3d & UNUSED(probe)) const
      __override {
    // WARNING This is totally arbitrary!
    return this->profile_.detEps;
  }
  virtual double coefficientCoulombDerivative_impl(
      const Eigen::Vector3d & UNUSED(direction), const Eigen::Vector3d & UNUSED(p1),
      const Eigen::Vector3d & UNUSED(p2)) const __override {
    return 0.0;
  }
  virtual double CoulombDerivative_impl(
      const Eigen::Vector3d & direction, const Eigen::Vector3d & p1,
      const Eigen::Vector3d & p2) const __override {
    return this->derivativeProbe(direction, p1, p2);
  }

  virtual double imagePotential_impl(const Eigen::Vector3d & UNUSED(source),
                                     const Eigen::Vector3d & UNUSED(probe)) const
      __override {
    return 0.0;
  }
  virtual double imagePotentialDerivative_impl(
      const Eigen::Vector3d & UNUSED(direction), const Eigen::Vector3d & UNUSED(p1),
      const Eigen::Vector3d & UNUSED(p2)) const __override {
    return 0.0;
  }

  virtual std::ostream & printObject(std::ostream & os) __override {
    os << "Green's function type: anisotropic liquid" << std::endl;
    os << this->profile_;
    return os;
  }
};

#endif // ANISOTROPICLIQUID_HPP
