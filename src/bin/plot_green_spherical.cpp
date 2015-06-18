#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "GeneralPurpose.hpp"
#include "LoggerInterface.hpp"
#include "OneLayerTanh.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"

std::vector<double> generate_grid(double, double, size_t);
void case1();
void case2();
void case3();
void case4();

int main()
{
    case1();
//    case2();
//    case3();
//    case4();

    return EXIT_SUCCESS;
}

std::vector<double> generate_grid(double min_val, double max_val, size_t nPoints)
{
    double step = (max_val - min_val) / nPoints;
    std::vector<double> grid;
    for (size_t i = 0; i < nPoints; ++i) {
        grid.push_back(min_val + i * step);
    }
    return grid;
}

void case1()
{
    size_t nPoints = 200;
    double epsInside = 80.0;
    double epsOutside = 2.0;
    Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
    double sphereRadius = 100.0;
    double width = 9.0;

    Eigen::Vector3d source = Eigen::Vector3d::Zero();
    Eigen::Vector3d probe = Eigen::Vector3d::Zero();

    double z_min = -20.0;
    double z_max = 20.0;
    std::vector<double> z = generate_grid(z_min, z_max, nPoints);
    double delta = 0.1;

    std::ofstream out;

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsInside, epsOutside, width, sphereRadius, sphereCenter);
    LOG(gf);

    out.open("gf_spherical_CASE1.dat");
    out.precision(16);
    for (size_t i = 0; i < nPoints; ++i) {
	source << 0.0, 0.0, z[i];
        for (size_t j = 0; j < nPoints; ++j) {
            probe << 0.0, 0.0, (z[j] + delta);
            out << '\t' << i
		<< '\t' << j
		<< '\t' << z[i]
                << '\t' << z[j]
                << '\t' << gf.Coulomb(source, probe)
                << '\t' << gf.imagePotential(source, probe)
                << '\t' << gf.kernelS(source, probe)
		<< '\t' << (1.0 / gf.coefficientCoulomb(source, probe)) << std::endl;
        }
    }
    out.close();
    LOG_TIME;

    /*
    UniformDielectric<AD_directional, CollocationIntegrator> gf_inside(epsInside);
    out.open("gf_uniform_inside_CASE1.dat");
    out << "#" << '\t' << "x_2" << '\t' << "z_2" << '\t';
    out << "gf_value" << '\t' << "derivativeProbe" << std::endl;
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            probe << x[i], 0.0, z[j];
            out << '\t' << x[i]
                << '\t' << z[j]
                << '\t' << gf_inside.kernelS(source, probe)
                << '\t' << gf_inside.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
        }
    }
    out.precision(16);
    out.close();

    UniformDielectric<AD_directional, CollocationIntegrator> gf_outside(epsOutside);
    out.open("gf_uniform_outside_CASE1.dat");
    out << "#" << '\t' << "x_2" << '\t' << "z_2" << '\t';
    out << "gf_value" << '\t' << "derivativeProbe" << std::endl;
    out.precision(16);
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            probe << x[i], 0.0, z[j];
            out << '\t' << x[i]
                << '\t' << z[j]
                << '\t' << gf_outside.kernelS(source, probe)
                << '\t' << gf_outside.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
        }
    }
    out.close();
    */
}
