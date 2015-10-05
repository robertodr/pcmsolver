/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "WaveletCavity.hpp"
#include "TestingMolecules.hpp"

/*! \class PWCSolver
 *  \test \b NH3 tests PWCSolver using ammonia and a wavelet cavity
 */
TEST_CASE("Piecewise constant wavelet solver for an ammonia molecule", "[solver][pwc][pwc_NH3]")
{
    LOG_INIT();
    Molecule molec = NH3();
    double probeRadius = 1.385; // Probe Radius for water
    int patchLevel = 4;
    double coarsity = 0.5;
    WaveletCavity cavity(spheres, probeRadius, patchLevel, coarsity);
    cavity.readCavity("molec_dyadic.dat");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator> gfInside;
    UniformDielectric<AD_directional, CollocationIntegrator> gfOutside(permittivity);
    int firstKind = 0;
    Compression comp(2.5, 2.5, 0.001);
    PWCSolver solver(comp, firstKind);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);
    cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_());

    double Ncharge = 7.0;
    double Hcharge = 1.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = computeMEP(molec, cavity.elements());
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    fake_asc = solver.computeCharge(fake_mep);
    double totalASC = - (Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    CAPTURE(totalASC - totalFakeASC);
    REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
}
