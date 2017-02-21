#ifndef TDPCMSOLVER_HPP
#define TDPCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class IGreensFunction;
class BoundaryIntegralOperator;

#include "Debye.hpp"

/*! \file TDPCMSolver.hpp
 *  \class TDPCMSolver
 *  \brief Abstract Base Class for time-dependent solvers inheritance hierarchy.
 *  \author Roberto Di Remigio
 *  \date 2015, 2016
 *  We use the Non-Virtual Interface idiom.
 */

class TDPCMSolver {
public:
  TDPCMSolver() {}
  /*! \brief Construct solver from two Green's functions
   *  \param[in] es static permittivity
   *  \param[in] ed dynamic permittivity
   *  \param[in] t  relaxation time
   */
  TDPCMSolver(double es, double ed, double t);
  virtual ~TDPCMSolver() {}

  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void buildSystemMatrix(const Cavity & cavity,
                         const IGreensFunction & gf_i,
                         const BoundaryIntegralOperator & op) {
    buildSystemMatrix_impl(cavity, gf_i, op);
  }
  /*! \brief Returns the ASC at time (t + dt) using a simple Euler integrator
   *  \param[in] dt the time step for the Euler integrator
   *  \param[in] MEP_current the vector containing the MEP at cavity points, at time
   * (t + dt)
   *  \param[in] MEP_previous the vector containing the MEP at cavity points, at time
   * t
   *  \param[in] ASC_previous the vector containing the ASC at cavity points, at time
   * t
   */
  Eigen::VectorXd propagateASC(double dt,
                               const Eigen::VectorXd & MEP_current,
                               const Eigen::VectorXd & MEP_previous,
                               const Eigen::VectorXd & ASC_previous) const {
    if (!built_)
      PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
    return propagateASC_impl(dt, MEP_current, MEP_previous, ASC_previous);
  }
  /*! \brief Returns the ASC at initial time
   *  \param[in] MEP the vector containing the MEP at cavity points at initial time
   *  Uses dynamic PCM matrix to compute the initial value for the ASC
   */
  Eigen::VectorXd initialValueASC(const Eigen::VectorXd & MEP) const {
    if (!built_)
      PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
    return initialValueASC_impl(MEP);
  }
  std::string printEnvironment();

  friend std::ostream & operator<<(std::ostream & os, TDPCMSolver & solver) {
    return solver.printSolver(os);
  }

protected:
  /*! Time-dependent isotropic, uniform Debye permittivity profile */
  Debye permittivity_;
  /*! Whether the system matrices have been built */
  bool built_;
  /*! Dynamic matrix */
  Eigen::MatrixXd A_;
  /*! Static matrix scaled by relaxation times */
  Eigen::MatrixXd B_;

  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  virtual void buildSystemMatrix_impl(const Cavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const BoundaryIntegralOperator & op) = 0;
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
      const Eigen::VectorXd & ASC_previous) const = 0;
  /*! \brief Returns the ASC at initial time
   *  \param[in] MEP the vector containing the MEP at cavity points at initial time
   *  Uses dynamic PCM matrix to compute the initial value for the ASC
   */
  virtual Eigen::VectorXd initialValueASC_impl(
      const Eigen::VectorXd & MEP) const = 0;
  virtual std::ostream & printSolver(std::ostream & os) = 0;
};

#endif // TDPCMSOLVER_HPP
