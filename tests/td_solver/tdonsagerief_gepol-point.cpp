#include "catch.hpp"

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>

#include "Config.hpp"

#include <Eigen/Core>

#include "utils/cnpy.hpp"
#include "bi_operators/Collocation.hpp"
#include "green/DerivativeTypes.hpp"
#include "cavity/GePolCavity.hpp"
#include "td_solver/TDSingleIEFSolver.hpp"
#include "TestingMolecules.hpp"
#include "green/Vacuum.hpp"
#include "green/UniformDielectric.hpp"
#include "solver/IEFSolver.hpp"
#include "PhysicalConstants.hpp"

using integrator::Collocation;

/*! \class TDOnsagerIEFSolver
 *  \test \b pointDipoleGePol tests TDOnsagerIEFSolver using a point dipole with a
 * GePol cavity
 */
TEST_CASE("Test solver for the time-dependent IEFPCM with the Onsager relaxation "
          "time, a point dipole and a GePol cavity",
          "[td_solver][tdief][tdonsagerief_gepol-point]") {
  // Set up cavity
  double radius = 1.181 * 1.10 / bohrToAngstrom();
  Molecule point = dummy<0>(radius);
  double area = 0.4;
  double probeRadius = 0.0;
  double minRadius = 100.0;
  GePolCavity cavity(point, area, probeRadius, minRadius);
  size_t size = cavity.size();

  double e_0 = 35.69;
  double e_d = 1.807;

  Vacuum<> gfInside;
  UniformDielectric<> gfOutside(e_d);
  Collocation biop;

  double tau = 2000; // 48.38 fs
  double dt = 0.2;   // 4.838 as
  double secondsToAU = 2.418884326509e-17;
  double total_time = 100e-15 / secondsToAU; // Total simulation time: 100 fs
  int steps = int(total_time / dt) + 1;      // Number of steps

  double tauOnsager = tau * (2 * e_d + 1) / (2 * e_0 + 1);
  TDSingleIEFSolver solver(e_0, e_d, tau, tauOnsager);
  solver.buildSystemMatrix(cavity, gfInside, biop);

  // Static solver
  bool symm = true;
  IEFSolver static_solver(symm);
  static_solver.buildSystemMatrix(cavity, gfInside, gfOutside, biop);

  // The point-like dipole is at the origin, this is the direction
  Eigen::Vector3d dipole = Eigen::Vector3d::UnitZ();
  Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
  for (size_t i = 0; i < size; ++i) {
    Eigen::Vector3d center = cavity.elementCenter(i);
    double distance = center.norm();
    fake_mep(i) = center.dot(dipole) / std::pow(distance, 3);
  }
  for (size_t i = 0; i < size; ++i) {
    INFO("fake_mep(" << i << ") = " << fake_mep(i));
  }
  Eigen::VectorXd fake_asc = static_solver.computeCharge(fake_mep);

  Eigen::VectorXd ASC_previous = solver.initialValueASC(fake_mep);
  for (size_t i = 0; i < size; ++i) {
    INFO("ASC_previous(" << i << ") = " << ASC_previous(i));
  }
  // Internal consistency check
  for (size_t i = 0; i < size; ++i) {
    REQUIRE(fake_asc(i) == Approx(ASC_previous(i)));
  }

  // Time-dependent CPCM reaction field
  std::vector<double> TDOnsagerIEF;
  TDOnsagerIEF.reserve(steps);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < steps; ++i) {
    TDOnsagerIEF.push_back(-fake_mep.dot(ASC_previous));
    asc = solver.propagateASC(dt, fake_mep, fake_mep, ASC_previous);
    ASC_previous = asc;
  }
  cnpy::NpyArray raw_TDOnsagerIEF_ref =
      cnpy::npy_load("tdonsagerief_collocation.npy");
  double * TDOnsagerIEF_ref = reinterpret_cast<double *>(raw_TDOnsagerIEF_ref.data);
  for (int i = 0; i < steps; ++i) {
    REQUIRE(TDOnsagerIEF[i] == Approx(TDOnsagerIEF_ref[i]));
  }
}
