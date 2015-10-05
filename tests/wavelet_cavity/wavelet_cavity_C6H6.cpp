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

TEST_CASE("Wavelet cavity for a benzene", "[wavelet_cavity][wavelet_cavity_C6H6]")
{
    LOG_INIT();
    Eigen::Vector3d C1( 5.274,  1.999, -8.568);
    Eigen::Vector3d C2( 6.627,  2.018, -8.209);
    Eigen::Vector3d C3( 7.366,  0.829, -8.202);
    Eigen::Vector3d C4( 6.752, -0.379, -8.554);
    Eigen::Vector3d C5( 5.399, -0.398, -8.912);
    Eigen::Vector3d C6( 4.660,  0.791, -8.919);
    Eigen::Vector3d H1( 4.704,  2.916, -8.573);
    Eigen::Vector3d H2( 7.101,  2.950, -7.938);
    Eigen::Vector3d H3( 8.410,  0.844, -7.926);
    Eigen::Vector3d H4( 7.322, -1.296, -8.548);
    Eigen::Vector3d H5( 4.925, -1.330, -9.183);
    Eigen::Vector3d H6( 3.616,  0.776, -9.196);

    std::vector<Sphere> spheres;
    Sphere sph1(C1, 1.53);
    Sphere sph2(C2, 1.53);
    Sphere sph3(C3, 1.53);
    Sphere sph4(C4, 1.53);
    Sphere sph5(C5, 1.53);
    Sphere sph6(C6, 1.53);

    Sphere sph7(H1, 1.06);
    Sphere sph8(H2, 1.06);
    Sphere sph9(H3, 1.06);
    Sphere sph10(H4, 1.06);
    Sphere sph11(H5, 1.06);
    Sphere sph12(H6, 1.06);

    spheres.push_back(sph1);
    spheres.push_back(sph2);
    spheres.push_back(sph3);
    spheres.push_back(sph4);
    spheres.push_back(sph5);
    spheres.push_back(sph6);
    spheres.push_back(sph7);
    spheres.push_back(sph8);
    spheres.push_back(sph9);
    spheres.push_back(sph10);
    spheres.push_back(sph11);
    spheres.push_back(sph12);

    double probeRadius = 1.385; // Probe Radius for water
    int patchLevel = 2;
    double coarsity = 0.5;
    cavity = WaveletCavity(spheres, probeRadius, patchLevel, coarsity);
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
     *  \test \b WaveletCavityC6H6Test_size tests Wavelet cavity size for benzene
     */
    SECTION("Test size")
    {
        size_t size = 6912;
        size_t actualSize = cavity.size();
        REQUIRE(size == actualSize);
    }

    /*! \class WaveletCavity
     *  \test \b WaveletCavityC6H6Test_area tests Wavelet cavity surface area for benzene
     */
    SECTION("Test surface area")
    {
        double area = 95.909894964414121;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
    }

    /*! \class WaveletCavity
     *  \test \b WaveletCavityC6H6Test_volume tests Wavelet cavity volume for benzene
     */
    SECTION("Test volume")
    {
        double volume = 69.622821450340595;
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
