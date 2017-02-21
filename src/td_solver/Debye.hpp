#ifndef DEBYE_HPP
#define DEBYE_HPP

#include <cmath>
#include <iosfwd>

#include "Config.hpp"

/*! \file Debye.hpp
 *  \struct Debye
 *  \brief A time-dependent Debye-type dielectric profile
 *  \author Roberto Di Remigio
 *  \date 2015
 */

struct Debye {
  Debye() {}
  Debye(double es, double ed, double t)
      : epsilonStatic(es), epsilonDynamic(ed), tau(t) {}
  /// Static dielectric constant
  double epsilonStatic;
  /// Dynamic dielectric constant
  double epsilonDynamic;
  /// Relaxation time
  double tau;
  /*! Return static Onsager factor
   *  \param[in] corr correction factor
   */
  double staticOnsager(double corr = 0.0) const {
    return -(epsilonStatic - 1) / (epsilonStatic + corr);
  }
  /*! Return dynamic Onsager factor
   *  \param[in] corr correction factor
   */
  double dynamicOnsager(double corr = 0.0) const {
    return -(epsilonDynamic - 1) / (epsilonDynamic + corr);
  }
  friend std::ostream & operator<<(std::ostream & os, Debye & obj) {
    os << "Profile functional form: Debye" << std::endl;
    os << "Static permittivity  = " << obj.epsilonStatic << std::endl;
    os << "Dynamic permittivity = " << obj.epsilonDynamic << std::endl;
    os << "Relaxation time      = " << obj.tau * AUToFemtoseconds() << " fs";
    return os;
  }
};

#endif // DEBYE_HPP
