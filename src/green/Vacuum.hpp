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

#ifndef VACUUM_HPP
#define VACUUM_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"

#include "GreenData.hpp"
#include "utils/ForId.hpp"
#include "utils/Factory.hpp"

/*! \file Vacuum.hpp
 *  \class Vacuum
 *  \brief Green's function for vacuum.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */
// TODO: * I don't think the ProfilePolicy parameter is needed in this case!
//       * can we use enable_if (or similar tricks) to avoid implementing useless
//       functions?

template <typename DerivativeTraits = AD_directional>
class Vacuum __final : public GreensFunction<DerivativeTraits, Uniform> {
public:
  Vacuum() : GreensFunction<DerivativeTraits, Uniform>() {
    this->profile_ = Uniform(1.0);
  }
  Vacuum(double fac) : GreensFunction<DerivativeTraits, Uniform>(fac) {
    this->profile_ = Uniform(1.0);
  }
  virtual ~Vacuum() {}

  friend std::ostream & operator<<(std::ostream & os, Vacuum & gf) {
    return gf.printObject(os);
  }

private:
  virtual DerivativeTraits operator()(DerivativeTraits * sp,
                                      DerivativeTraits * pp) const __override {
    return (1 / distance(sp, pp));
  }
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const __override {
    return this->derivativeProbe(direction, p1, p2);
  }

  virtual KernelS exportKernelS_impl() const __override {
    return pcm::bind(&Vacuum<DerivativeTraits>::kernelS, *this, pcm::_1, pcm::_2);
  }
  virtual KernelD exportKernelD_impl() const __override {
    return pcm::bind(&Vacuum<DerivativeTraits>::kernelD, *this, pcm::_1, pcm::_2,
                     pcm::_3);
  }

  virtual double singleLayer_impl(const Element & e) const __override {
    return integrator::diagonalSi(e.area(), this->factor_);
  }
  virtual double doubleLayer_impl(const Element & e) const __override {
    return integrator::diagonalDi(e.area(), e.sphere().radius, this->factor_);
  }

  virtual std::ostream & printObject(std::ostream & os) __override {
    os << "Green's function type: vacuum";
    return os;
  }
};

namespace {
struct buildVacuum {
  template <typename T> IGreensFunction * operator()(const greenData & data) {
    return new Vacuum<T>(data.scaling);
  }
};

IGreensFunction * createVacuum(const greenData & data) {
  buildVacuum build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}
const std::string VACUUM("VACUUM");
const bool registeredVacuum =
    Factory<IGreensFunction, greenData>::TheFactory().registerObject(VACUUM,
                                                                     createVacuum);
}

#endif // VACUUM_HPP
