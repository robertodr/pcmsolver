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

#define BOOST_TEST_MODULE IEFSolverpointChargeGePolSymmetry

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/filesystem.hpp>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "Symmetry.hpp"

namespace fs = boost::filesystem;

/*! \class IEFSolver
 *  \test \b pointChargeGePolC1 tests IEFSolver using a point charge with a GePol cavity in C1 symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC1)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // C1
    Symmetry group = buildGroup(0, 0, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c1");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolCs tests IEFSolver using a point charge with a GePol cavity in Cs symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolCs)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // Cs as generated by Oyz
    Symmetry group = buildGroup(1, 1, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.cs");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolC2 tests IEFSolver using a point charge with a GePol cavity in C2 symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC2)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // C2 as generated by C2z
    Symmetry group = buildGroup(1, 3, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolCi tests IEFSolver using a point charge with a GePol cavity in Ci symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolCi)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // Ci as generated by i
    Symmetry group = buildGroup(1, 7, 0, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.ci");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolC2h tests IEFSolver using a point charge with a GePol cavity in C2h symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC2h)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // C2h as generated by Oxy and i
    Symmetry group = buildGroup(2, 4, 7, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2h");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolD2 tests IEFSolver using a point charge with a GePol cavity in D2 symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolD2)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // D2 as generated by C2z and C2x
    Symmetry group = buildGroup(2, 3, 6, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.d2");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolC2v tests IEFSolver using a point charge with a GePol cavity in C2v symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC2v)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // C2v as generated by Oyz and Oxz
    Symmetry group = buildGroup(2, 1, 2, 0);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class IEFSolver
 *  \test \b pointChargeGePolD2h tests IEFSolver using a point charge with a GePol cavity in D2h symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolD2h)
{
    // Set up cavity
    Eigen::Vector3d N(0.0, 0.0, 0.0);
    std::vector<Sphere> spheres;
    Sphere sph1(N, 2.929075493);
    spheres.push_back(sph1);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    // D2h as generated by Oxy, Oxz and Oyz
    Symmetry group = buildGroup(3, 4, 2, 1);
    GePolCavity cavity(spheres, area, probeRadius, minRadius, group);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.d2h");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}
