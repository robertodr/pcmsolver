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

#ifndef UNIFORMDIELECTRIC_HPP
#define UNIFORMDIELECTRIC_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "GreensFunction.hpp"
#include "dielectric_profile/Uniform.hpp"

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric
 *  \brief Green's function for uniform dielectric.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */

template <typename DerivativeTraits = AD_directional>
class UniformDielectric __final : public GreensFunction<DerivativeTraits, Uniform> {
public:
  UniformDielectric(double eps) : GreensFunction<DerivativeTraits, Uniform> {
    this->profile_ = Uniform(eps);
  }
  virtual ~UniformDielectric() {}

  double epsilon() const { return this->profile_.epsilon; }

  friend std::ostream & operator<<(std::ostream & os, UniformDielectric & gf) {
    return gf.printObject(os);
  }

private:
  virtual DerivativeTraits operator()(DerivativeTraits * sp,
                                      DerivativeTraits * pp) const __override {
    return 1 / (this->profile_.epsilon * distance(sp, pp));
  }
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const __override {
    return this->profile_.epsilon * (this->derivativeProbe(direction, p1, p2));
  }

  virtual KernelS exportKernelS_impl() const __override {
    return pcm::bind(&UniformDielectric<DerivativeTraits, IntegratorPolicy>::kernelS,
                     *this, pcm::_1, pcm::_2);
  }
  virtual KernelD exportKernelD_impl() const __override {
    return pcm::bind(&UniformDielectric<DerivativeTraits, IntegratorPolicy>::kernelD,
                     *this, pcm::_1, pcm::_2, pcm::_3);
  }

  virtual double coefficient_impl(const Eigen::Vector3d & UNUSED(source),
                                  const Eigen::Vector3d & UNUSED(probe)) const
      __override {
    return this->profile_.epsilon;
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
    os << "Green's function type: uniform dielectric" << std::endl;
    os << this->profile_;
    return os;
  }
};

#endif // UNIFORMDIELECTRIC_HPP
