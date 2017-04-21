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

#include "ddPCM.hpp"

namespace pcm {
ddPCM::ddPCM(const Molecule & m) {
  int size = m.spheres().size();
  double * xs = new double[size];
  double * ys = new double[size];
  double * zs = new double[size];
  double * rs = new double[size];
  for (int i = 0; i < size; ++i) {
    xs(i) = m.spheres(i).center(0);
    ys(i) = m.spheres(i).center(1);
    zs(i) = m.spheres(i).center(2);
    rs(i) = m.spheres(i).radius;
  }
  ddinit(size, xs, ys, zs, rs);
  delete[] xs, ys, zs, rs;
}

ddPCM::~ddPCM() { memfree(); }
}
