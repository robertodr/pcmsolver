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

#ifndef TDSOLVER_HPP
#define TDSOLVER_HPP

#include "Config.hpp"

#include "ITDSolver.hpp"
#include "TDCPCMSolver.hpp"
#include "TDIEFSolver.hpp"
#include "TDSingleIEFSolver.hpp"
#include "utils/Factory.hpp"

/*!
 * \file TDSolver.hpp
 * \brief Top-level include file for time-dependent solvers
 * \author Roberto Di Remigio
 * \date 2017
 *
 * Includes all time-dependent solver-related headers and defines the bootstrap
 *function
 * for the Factory<ITDSolver, TDSolverData>
 */

namespace pcm {
namespace td_solver {
namespace detail {
typedef pcm::function<ITDSolver *(const TDSolverData &)> CreateTDSolver;
} // namespace detail

inline Factory<detail::CreateTDSolver> bootstrapFactory() {
  Factory<detail::CreateTDSolver> factory_;

  factory_.subscribe("TDCPCM", createTDCPCMSolver);
  factory_.subscribe("TDIEF", createTDIEFSolver);
  factory_.subscribe("TDSINGLEIEF", createTDSingleIEFSolver);
  factory_.subscribe("TDONSAGERIEF", createTDOnsagerIEFSolver);

  return factory_;
}
} // namespace td_solver
} // namespace pcm

#endif // TDSOLVER_HPP
