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

#define BOOST_TEST_MODULE GePolCavityC2vTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <cstdio>

#include "Config.hpp"

#include <Eigen/Core>


#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"
#include "TestingMolecules.hpp"


// Test C2v symmetry with addition of extra spheres enabled
struct GePolCavityC2vAddTest {
    GePolCavity cavity;
    GePolCavityC2vAddTest() { SetUp(); }
    void SetUp() {
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
	Molecule molec = H3<5>();
        cavity = GePolCavity(molec, area, probeRadius, minRadius);
        cavity.saveCavity("h3+_c2v.npz");
        std::rename("PEDRA.OUT", "PEDRA.OUT.c2v");
        std::rename("cavity.off", "cavity.off.c2v");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityC2vAdd, GePolCavityC2vAddTest)

/*! \class GePolCavity
 *  \test \b GePolCavityC2vAddTest_size tests GePol cavity size for H3+ in C2v symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC2vAddTest)
{
    int size = 312;
    size_t actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vAddTest_irreducible_size tests GePol cavity irreducible size for H3+ in C2v symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityC2vAddTest)
{
    int size = 78;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vAddTest_area tests GePol cavity surface area for H3+ in C2v symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC2vAddTest)
{
    double area = 178.74700256128352;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vAddTest_volume tests GePol cavity volume for H3+ in C2v symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC2vAddTest)
{
    double volume = 196.47360294559090;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( size_t i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()

// Test C2v symmetry without addition of extra spheres enabled
struct GePolCavityC2vTest {
    GePolCavity cavity;
    GePolCavityC2vTest() { SetUp(); }
    void SetUp() {
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 100.0 / convertBohrToAngstrom;
	Molecule molec = H3<5>();
        cavity = GePolCavity(molec, area, probeRadius, minRadius);
        cavity.saveCavity("h3+_c2v_noadd.npz");
        std::rename("PEDRA.OUT", "PEDRA.OUT.c2v_noadd");
        std::rename("cavity.off", "cavity.off.c2v_noadd");
    }
};

BOOST_FIXTURE_TEST_SUITE(GePolCavityC2v, GePolCavityC2vTest)

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_size tests GePol cavity size for H3+ in C2v symmetry without added spheres
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC2vTest)
{
    int size = 288;
    size_t actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_irreducible_size tests GePol cavity irreducible size for H3+ in C2v symmetry without added spheres
 */
BOOST_FIXTURE_TEST_CASE(irreducible_size, GePolCavityC2vTest)
{
    int size = 72;
    int actualSize = cavity.irreducible_size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_area tests GePol cavity surface area for H3+ in C2v symmetry without added spheres
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC2vTest)
{
    double area = 181.87043332808548;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC2vTest_volume tests GePol cavity volume for H3+ in C2v symmetry without added spheres
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC2vTest)
{
    double volume = 192.48281460140359;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( size_t i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END()
