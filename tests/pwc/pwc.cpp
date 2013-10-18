#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "PWCSolver.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "WaveletCavity.hpp"
#include "gtest/gtest.h"

TEST(PWCSolver, pointCharge) 
{
	// Set up cavity
	Eigen::Vector3d N(0.0, 0.0, 0.0); 	
	std::vector<Sphere> spheres;      	
	Sphere sph1(N, 2.929075493);      	
	spheres.push_back(sph1);
	double probeRadius = 1.385;
	int patchLevel = 2;
	double coarsity = 0.5;
	WaveletCavity cavity(spheres, probeRadius, patchLevel, coarsity);
	cavity.readCavity("molec_dyadic.dat");
	// The point charge is located at the origin.
	// The potential at cavity point s_I is Q/|s_I|
	double permittivity = 78.39;
	Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
	UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
	int firstKind = 0;
	PWCSolver solver(gfInside, gfOutside, firstKind);
	solver.buildSystemMatrix(cavity);
	cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);

	double charge = 8.0;
	int size = cavity.size();
	Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		Eigen::Vector3d center = cavity.getElementCenter(i);
		double distance = center.norm();
		fake_mep(i) = charge / distance; 
	}
	// The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
	Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
	solver.compCharge(fake_mep, fake_asc);
	double totalASC = - charge * (permittivity - 1) / permittivity;
	double totalFakeASC = fake_asc.sum();
	std::cout << "totalASC -totalFakeASC = " << totalASC - totalFakeASC << std::endl;
	EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}