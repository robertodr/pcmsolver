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

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "WaveletCavity.hpp"
#include "Molecule.hpp"
#include "PhysicalConstants.hpp"
#include "TestingMolecules.hpp"

TEST_CASE("Wavelet cavity for an ammonia molecule", "[wavelet_cavity][wavelet_cavity_NH3]")
{
    LOG_INIT();
    Molecule molec = NH3();
    double probeRadius = 1.385; // Probe Radius for water
    int patchLevel = 2;
    double coarsity = 0.5;
    WaveletCavity cavity = WaveletCavity(molec.spheres(), probeRadius, patchLevel, coarsity);
    cavity.readCavity("molec_dyadic.dat");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator> gfInside;
    UniformDielectric<AD_directional, CollocationIntegrator> gfOutside(permittivity);
    int firstKind = 0;
    Compression comp(2.5, 2.5, 0.001);
    PWCSolver solver(comp, firstKind);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);
    cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_());

    /*! \class WaveletCavity
     *  \test \b WaveletCavityNH3Test_size tests Wavelet cavity size for ammonia
     */
    SECTION("Test size")
    {
        size_t size = 4288;
        size_t actualSize = cavity.size();
        REQUIRE(size == actualSize);
    }

    /*! \class WaveletCavity
     *  \test \b WaveletCavityNH3Test_area tests Wavelet cavity surface area for ammonia
     */
    SECTION("Test surface area")
    {
        double area = 146.41490284471513;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
    }

    /*! \class WaveletCavity
     *  \test \b WaveletCavityNH3Test_volume tests Wavelet cavity volume for ammonia
     */
    SECTION("Test volume")
    {
        double volume = 153.3346491517182;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i ) {
            actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                        i));
        }
        actualVolume /= 3;
        REQUIRE(volume == Approx(actualVolume));
    }
}
