/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#ifndef ITDSOLVER_HPP
#define ITDSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
class ICavity;
class IGreensFunction;
class IBoundaryIntegralOperator;
} // namespace pcm

#include "Debye.hpp"

/*! \file ITDSolver.hpp
 *  \class ITDSolver
 *  \brief Abstract Base Class for time-dependent solvers inheritance hierarchy.
 *  \author Roberto Di Remigio
 *  \date 2017
 *  We use the Non-Virtual Interface idiom.
 */

namespace pcm {
using td_solver::Debye;
class ITDSolver {
public:
  ITDSolver() {}
  /*! \brief Construct solver from two Green's functions
   *  \param[in] es static permittivity
   *  \param[in] ed dynamic permittivity
   *  \param[in] t  relaxation time
   */
  ITDSolver(double es, double ed, double t);
  virtual ~ITDSolver() {}

  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void buildSystemMatrix(const ICavity & cavity,
                         const IGreensFunction & gf_i,
                         const IBoundaryIntegralOperator & op) {
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
      PCMSOLVER_ERROR("PCM matrix not calculated yet");
    return propagateASC_impl(dt, MEP_current, MEP_previous, ASC_previous);
  }
  /*! \brief Returns the ASC at initial time
   *  \param[in] MEP the vector containing the MEP at cavity points at initial time
   *  Uses dynamic PCM matrix to compute the initial value for the ASC
   */
  Eigen::VectorXd initialValueASC(const Eigen::VectorXd & MEP) const {
    if (!built_)
      PCMSOLVER_ERROR("PCM matrix not calculated yet");
    return initialValueASC_impl(MEP);
  }
  std::string printEnvironment();

  friend std::ostream & operator<<(std::ostream & os, ITDSolver & solver) {
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
  virtual void buildSystemMatrix_impl(const ICavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const IBoundaryIntegralOperator & op) = 0;
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
} // namespace pcm

#endif // ITDSOLVER_HPP
