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

#ifndef TDSINGLEIEFSOLVER_HPP
#define TDSINGLEIEFSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class IGreensFunction;
class BoundaryIntegralOperator;

#include "TDPCMSolver.hpp"

/*! \file TDSingleIEFSolver.hpp
 *  \class TDSingleIEFSolver
 *  \brief Time-dependent solver for isotropic IEF with a single relaxation time
 *  \author Roberto Di Remigio
 *  \date 2015, 2016
 */

class TDSingleIEFSolver : public TDPCMSolver {
public:
  TDSingleIEFSolver() {}
  /*! \brief Construct solver from two Green's functions
   *  \param[in] es static permittivity
   *  \param[in] ed dynamic permittivity
   *  \param[in] t  relaxation time
   *  \param[in] tau IEF relaxation time
   */
  TDSingleIEFSolver(double es, double ed, double t, double tau);
  virtual ~TDSingleIEFSolver() {}
  friend std::ostream & operator<<(std::ostream & os, TDSingleIEFSolver & solver) {
    return solver.printSolver(os);
  }

protected:
  /*! IEF relaxation time */
  double tauIEF_;

private:
  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
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

#endif // TDSINGLEIEFSOLVER_HPP
