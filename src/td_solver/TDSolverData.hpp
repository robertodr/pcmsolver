/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

#pragma once

#include "Config.hpp"

namespace pcm {
/*! \struct TDSolverData
 *  \brief Contains all data defined from user input in the solver section.
 */
struct TDSolverData {
  /*! Static permittivity */
  double epsilonStatic;
  /*! Dynamic permittivity */
  double epsilonDynamic;
  /*! Relaxatione time */
  double tau;
  /*! The correction factor to be use in a CPCM calculation */
  double correction;
  /*! IEF relaxation time */
  double tauIEF;
  /*! Whether to use Cholesky decomposition */
  bool cholesky;
  /*! Whether to initialize time-evolution with dynamic ASC */
  bool initWithDynamic;
  /*! Whether the structure was initialized with user input or not */
  bool empty;

  TDSolverData() { empty = true; }
  TDSolverData(double es,
               double ed,
               double t,
               double corr,
               double tau,
               bool dyn,
               bool chol)
      : epsilonStatic(es),
        epsilonDynamic(ed),
        tau(t),
        correction(corr),
        tauIEF(tau),
        cholesky(chol),
        initWithDynamic(dyn) {
    empty = false;
  }
};
} // namespace pcm
