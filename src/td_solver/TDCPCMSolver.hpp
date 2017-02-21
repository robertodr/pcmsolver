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

#ifndef TDCPCMSOLVER_HPP
#define TDCPCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
class ICavity;
class IGreensFunction;
class IBoundaryIntegralOperator;
struct TDSolverData;
} // namespace pcm

#include "ITDSolver.hpp"

/*! \file TDCPCMSolver.hpp
 *  \class TDCPCMSolver
 *  \brief Time-dependent solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2015-2017
 */

namespace pcm {
namespace td_solver {
class TDCPCMSolver : public ITDSolver {
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
  virtual void buildSystemMatrix_impl(const ICavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const IBoundaryIntegralOperator & op)
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

ITDSolver * createTDCPCMSolver(const TDSolverData & data);
} // namespace td_solver
} // namespace pcm

#endif // TDCPCMSOLVER_HPP
