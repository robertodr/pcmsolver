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

#ifndef TDIEFSOLVER_HPP
#define TDIEFSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class IGreensFunction;
class BoundaryIntegralOperator;

#include "TDPCMSolver.hpp"

/*! \file TDIEFSolver.hpp
 *  \class TDIEFSolver
 *  \brief Time-dependent solver for the isotropic IEF model
 *  \author Roberto Di Remigio
 *  \date 2015
 */

class TDIEFSolver : public TDPCMSolver {
public:
  TDIEFSolver() {}
  /*! \brief Construct time-dependent solver
   *  \param[in] es static permittivity
   *  \param[in] ed dynamic permittivity
   *  \param[in] t  relaxation time
   *  \param[in] cholesky whether Cholesky or Lowdin decomposition is used
   */
  TDIEFSolver(double es, double ed, double t, bool cholesky);
  virtual ~TDIEFSolver() {}
  /*! Return i-th element of the static K matrix
   *  \param[in] i index
   */
  double K_0(int i) const { return K_0_(i); }
  /*! Return i-th element of the dynamic K matrix
   *  \param[in] i index
   */
  double K_d(int i) const { return K_d_(i); }
  /*! Return i-th element of the Lambda matrix
   *  \param[in] i index
   */
  double Lambda(int i) const { return Lambda_(i); }
  /*! Return i-th element of the relaxation times matrix
   *  \param[in] i index
   *  \note Relaxation times are expressed in atomic units
   */
  double tau(int i) const { return tau_(i); }
  /*! Return static K matrix */
  const Eigen::VectorXd & K_0() const { return K_0_; }
  /*! Return dynamic K matrix */
  const Eigen::VectorXd & K_d() const { return K_d_; }
  /*! Return Lambda matrix */
  const Eigen::VectorXd & Lambda() const { return Lambda_; }
  /*! Return relaxation times matrix
   *  \note Relaxation times are expressed in atomic units
   */
  const Eigen::VectorXd & tau() const { return tau_; }
  friend std::ostream & operator<<(std::ostream & os, TDIEFSolver & solver) {
    return solver.printSolver(os);
  }

private:
  /*! Whether Cholesky decomposition is used to solve the generalized eigenvalue
   * problem */
  bool useCholesky_;
  /*! Diagonalized PCM matrix, static */
  Eigen::VectorXd K_0_;
  /*! Diagonalized PCM matrix, dynamic */
  Eigen::VectorXd K_d_;
  /*! Diagonalized S^-1/2DAS^1/2 matrix */
  Eigen::VectorXd Lambda_;
  /*! Diagonal relaxation times matrix */
  Eigen::VectorXd tau_;
  /*! Relaxation times (used in the propagator) */
  Eigen::MatrixXd C_;
  /*! Dynamic isotropic IEF matrix (used for the initial value of the ASC) */
  Eigen::MatrixXd PCMMatrix_;

  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  virtual void buildSystemMatrix_impl(const Cavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const BoundaryIntegralOperator & op)
      __override;
  /*! \brief Calculation of the PCM matrix, using Lowdin symmetric orthogonalization
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void systemMatrix_Lowdin(const Cavity & cavity,
                           const IGreensFunction & gf_i,
                           const BoundaryIntegralOperator & op);
  /*! \brief Calculation of the PCM matrix, using Cholesky orthogonalization
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i   Green's function inside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void systemMatrix_Cholesky(const Cavity & cavity,
                             const IGreensFunction & gf_i,
                             const BoundaryIntegralOperator & op);
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

#endif // TDIEFSOLVER_HPP
