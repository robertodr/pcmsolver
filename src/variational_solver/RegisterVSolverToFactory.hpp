/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef REGISTER_VSOLVER_TO_FACTORY_HPP
#define REGISTER_VSOLVER_TO_FACTORY_HPP

#include <string>

#include "Factory.hpp"
#include "SolverData.hpp"
#include "VCPCMSolver.hpp"
#include "VIEFSolver.hpp"

/*! \file RegisterVSolverToFactory.hpp
 *  \brief Register each variational solver to the factory.
 *  \author Roberto Di Remigio
 *  \date 2015
 */

// TODO These static casts are ugly, think of a better solution?

namespace
{
    VPCMSolver * createVCPCMSolver(const solverData & data)
    {
        return new VCPCMSolver(static_cast<VPCMSolver::GuessType>(data.guess), data.correction);
    }
    const std::string VCPCMSOLVER("VCPCM");
    const bool registeredVCPCMSolver =
        Factory<VPCMSolver, solverData>::TheFactory().registerObject(VCPCMSOLVER, createVCPCMSolver);
}

namespace
{
    VPCMSolver * createVIEFSolver(const solverData & data)
    {
        return new VIEFSolver(static_cast<VPCMSolver::GuessType>(data.guess));
    }
    const std::string VIEFSOLVER("VIEFPCM");
    const bool registeredVIEFSolver =
        Factory<VPCMSolver, solverData>::TheFactory().registerObject(VIEFSOLVER, createVIEFSolver);
}

#endif // REGISTER_VSOLVER_TO_FACTORY_HPP
