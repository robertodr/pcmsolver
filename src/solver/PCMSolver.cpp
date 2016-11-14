/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators.
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

#include "PCMSolver.hpp"

#include "Config.hpp"

#include "green/IGreensFunction.hpp"
#include "green/dielectric_profile/ProfileTypes.hpp"

PCMSolver::EnvironmentType PCMSolver::environment(const IGreensFunction & gf_i,
                                                  const IGreensFunction & gf_o) {
  PCMSolver::EnvironmentType retval = PCMSolver::AnisotropicEnv;
  bool isotropic = (gf_i.uniform() && gf_o.uniform());
  // Determine if we are in the regular or flipped regime
  if (isotropic) {
    bool flipped = ((profiles::epsilon(gf_i.permittivity()) != 1.0) &&
                    (profiles::epsilon(gf_o.permittivity()) == 1.0));
    flipped ? retval = PCMSolver::FlippedIsotropic : retval =
                                                         PCMSolver::RegularIsotropic;
  }
  return retval;
}
