#ifndef TDCPCMSOLVER_HPP
#define TDCPCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class IGreensFunction;
class BoundaryIntegralOperator;

#include "TDPCMSolver.hpp"

/*! \file TDCPCMSolver.hpp
 *  \class TDCPCMSolver
 *  \brief Time-dependent solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2015
 */

class TDCPCMSolver : public TDPCMSolver {
public:
  TDCPCMSolver() {}
  /*! \brief Construct solver from two Green's functions
   *  \param[in] es static permittivity
   *  \param[in] ed dynamic permittivity
   *  \param[in] t  relaxation time
   *  \param[in] corr factor to correct the conductor results
   */
  TDCPCMSolver(double es, double ed, double t, double corr);
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
  virtual void buildSystemMatrix_impl(const Cavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const BoundaryIntegralOperator & op)
      __override;
  /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
   *  \param[in] dt the time step for the Euler integrator
   *  \param[in] MEP_current the vector containing the MEP at cavity points, at time
   * (t + dt)
   *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time
   * t
   *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time
   * t
   */
  virtual Eigen::VectorXd propagateASC_impl(
      double dt,
      const Eigen::VectorXd & MEP_current,
      const Eigen::VectorXd & MEP_previous,
      const Eigen::VectorXd & ASC_previous) const __override;
  /*! \brief Returns the ASC at initial time
   *  \param[in] MEP the vector containing the MEP at cavity points at initial time
   *  Uses dynamic PCM matrix to compute the initial value for the ASC
   */
  virtual Eigen::VectorXd initialValueASC_impl(const Eigen::VectorXd & MEP) const
      __override;
  virtual std::ostream & printSolver(std::ostream & os) __override;
};

#endif // TDCPCMSOLVER_HPP
