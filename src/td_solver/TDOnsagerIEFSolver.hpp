#ifndef TDONSAGERIEFSOLVER_HPP
#define TDONSAGERIEFSOLVER_HPP

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
#include "TDSingleIEFSolver.hpp"

/*! \file TDOnsagerIEFSolver.hpp
 *  \class TDOnsagerIEFSolver
 *  \brief Time-dependent solver for isotropic IEF with a single Onsager relaxation
 * time
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation
 * of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class TDOnsagerIEFSolver
    : public TDSingleIEFSolver<DerivativeTraits, IntegratorPolicy> {
public:
  TDOnsagerIEFSolver() {}
  /*! \brief Construct solver
   *  \param[in] es static permittivity
   *  \param[in] ed dynamic permittivity
   *  \param[in] t  relaxation time
   */
  TDOnsagerIEFSolver(double es, double ed, double t)
      : TDSingleIEFSolver<DerivativeTraits, IntegratorPolicy>(
            es,
            ed,
            t,
            (t * (2 * ed + 1) / (2 * es + 1))) {}
  virtual ~TDOnsagerIEFSolver() {}
  friend std::ostream & operator<<(std::ostream & os, TDOnsagerIEFSolver & solver) {
    return solver.printSolver(os);
  }

private:
  virtual std::ostream & printSolver(std::ostream & os) {
    os << "Solver Type: IEFPCM, isotropic" << std::endl;
    os << "Onsager relaxation time = " << this->tauIEF_;
    return os;
  }
};

#endif // TDONSAGERIEFSOLVER_HPP
