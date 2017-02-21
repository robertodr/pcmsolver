#include "TDSingleIEFSolver.hpp"

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "bi_operators/BoundaryIntegralOperator.hpp"
#include "cavity/Cavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "Debye.hpp"
#include "TDPCMSolver.hpp"
#include "TDSolverData.hpp"
#include "utils/Factory.hpp"

TDSingleIEFSolver::TDSingleIEFSolver(double es, double ed, double t, double tau)
    : TDPCMSolver(es, ed, t), tauIEF_(tau) {}

void TDSingleIEFSolver::buildSystemMatrix_impl(const Cavity & cavity,
                                               const IGreensFunction & gf_i,
                                               const BoundaryIntegralOperator & op) {
  // The total size of the cavity
  int cavitySize = cavity.size();
  // Identity matrix
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);
  double e_d = permittivity_.epsilonDynamic;
  double e_0 = permittivity_.epsilonStatic;
  // Compute S and D on the whole cavity, regardless of symmetry
  Eigen::MatrixXd S = op.computeS(cavity, gf_i);
  Eigen::MatrixXd D = op.computeD(cavity, gf_i);
  Eigen::MatrixXd A = cavity.elementArea().asDiagonal();
  // Form A_ (the dynamic matrix)
  double f_d = (e_d + 1.0) / (e_d - 1.0);
  A_ = 2 * M_PI * f_d * S - D * A * S;
  Eigen::FullPivLU<Eigen::MatrixXd> A__LU(A_);
  if (!(A__LU.isInvertible()))
    PCMSOLVER_ERROR("A_ matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  A_ = -A__LU.inverse();
  A_ *= (2 * M_PI * Id - D * A);
  hermitivitize(A_);
  // Form B_ (the static matrix)
  double f_0 = (e_0 + 1.0) / (e_0 - 1.0);
  B_ = 2 * M_PI * f_0 * S - D * A * S;
  Eigen::FullPivLU<Eigen::MatrixXd> B__LU(B_);
  if (!(B__LU.isInvertible()))
    PCMSOLVER_ERROR("B_ matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  B_ = -B__LU.inverse();
  B_ *= (2 * M_PI * Id - D * A);
  hermitivitize(B_);
  built_ = true;
}

Eigen::VectorXd TDSingleIEFSolver::propagateASC_impl(
    double dt,
    const Eigen::VectorXd & MEP_current,
    const Eigen::VectorXd & MEP_previous,
    const Eigen::VectorXd & ASC_previous) const {
  double factor = dt / tauIEF_;
  return (A_ * (MEP_current - MEP_previous) + factor * B_ * MEP_previous -
          factor * ASC_previous + ASC_previous);
}

Eigen::VectorXd TDSingleIEFSolver::initialValueASC_impl(
    const Eigen::VectorXd & MEP) const {
  return (A_ * MEP);
}

std::ostream & TDSingleIEFSolver::printSolver(std::ostream & os) {
  os << "Solver Type: IEFPCM, isotropic" << std::endl;
  os << "IEF relaxation time = " << tauIEF_;
  return os;
}

namespace {
TDPCMSolver * createTDSingleIEFSolver(const TDSolverData & data) {
  return new TDSingleIEFSolver(
      data.epsilonStatic, data.epsilonDynamic, data.tau, data.tauIEF);
}
const std::string TDSINGLEIEFSOLVER("TDSINGLEIEF");
const bool registeredTDSingleIEFSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDSINGLEIEFSOLVER,
        createTDSingleIEFSolver);
}

namespace {
TDPCMSolver * createTDOnsagerIEFSolver(const TDSolverData & data) {
  double tauOnsager =
      data.tau * (2 * data.epsilonDynamic + 1) / (2 * data.epsilonStatic + 1);
  return new TDSingleIEFSolver(
      data.epsilonStatic, data.epsilonDynamic, data.tau, tauOnsager);
}
const std::string TDONSAGERIEFSOLVER("TDONSAGERIEF");
const bool registeredTDOnsagerIEFSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDONSAGERIEFSOLVER,
        createTDOnsagerIEFSolver);
}
