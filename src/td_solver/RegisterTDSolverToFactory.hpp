#ifndef REGISTERTDSOLVERTOFACTORY_HPP
#define REGISTERTDSOLVERTOFACTORY_HPP

#include <string>

#include "Config.hpp"

#include "green/DerivativeTypes.hpp"
#include "utils/Factory.hpp"
#include "utils/ForId.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "bi_operators/IntegratorTypes.hpp"
#include "TDCPCMSolver.hpp"
#include "TDCPCMIterativeSolver.hpp"
#include "TDIEFSolver.hpp"
#include "TDOnsagerIEFSolver.hpp"
#include "TDSingleIEFSolver.hpp"
#include "TDSolverData.hpp"

/*! \file RegisterTDSolverToFactory.hpp
 *  \brief Register each TD solver to the factory.
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  This file collects all the creational functions needed for the creation
 *  of the solver objects by means of the factory method.
 *  Originally, each of them was in the same file as the respective class.
 *  This, however, lead to intricate inclusion dependencies.
 */

namespace {
struct buildTDCPCMSolver {
  template <typename T, typename U>
  TDPCMSolver * operator()(const TDSolverData & data) {
    return new TDCPCMSolver<T, U>(
        data.epsilonStatic, data.epsilonDynamic, data.tau, data.correction);
  }
};

TDPCMSolver * createTDCPCMSolver(const TDSolverData & data) {
  buildTDCPCMSolver build;
  return for_id<derivative_types, integrator_types, TDPCMSolver>(
      build, data, data.howDerivative, data.howIntegrator);
}
const std::string TDCPCMSOLVER("TDCPCM");
const bool registeredTDCPCMSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDCPCMSOLVER,
        createTDCPCMSolver);
}

namespace {
struct buildTDCPCMIterativeSolver {
  template <typename T, typename U>
  TDPCMIterativeSolver * operator()(const TDSolverData & data) {
    return new TDCPCMIterativeSolver<T, U>(
        data.epsilonStatic, data.epsilonDynamic, data.tau, data.correction);
  }
};

TDPCMIterativeSolver * createTDCPCMIterativeSolver(const TDSolverData & data) {
  buildTDCPCMIterativeSolver build;
  return for_id<derivative_types, integrator_types, TDPCMIterativeSolver>(
      build, data, data.howDerivative, data.howIntegrator);
}
const std::string TDCPCMIterativeSOLVER("TDCPCMITERATIVE");
const bool registeredTDCPCMIterativeSolver =
    Factory<TDPCMIterativeSolver, TDSolverData>::TheFactory().registerObject(
        TDCPCMIterativeSOLVER,
        createTDCPCMIterativeSolver);
}

namespace {
struct buildTDIEFSolver {
  template <typename T, typename U>
  TDPCMSolver * operator()(const TDSolverData & data) {
    return new TDIEFSolver<T, U>(
        data.epsilonStatic, data.epsilonDynamic, data.tau, data.cholesky);
  }
};

TDPCMSolver * createTDIEFSolver(const TDSolverData & data) {
  buildTDIEFSolver build;
  return for_id<derivative_types, integrator_types, TDPCMSolver>(
      build, data, data.howDerivative, data.howIntegrator);
}
const std::string TDIEFSOLVER("TDIEF");
const bool registeredTDIEFSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDIEFSOLVER,
        createTDIEFSolver);
}

namespace {
struct buildTDOnsagerIEFSolver {
  template <typename T, typename U>
  TDPCMSolver * operator()(const TDSolverData & data) {
    return new TDOnsagerIEFSolver<T, U>(
        data.epsilonStatic, data.epsilonDynamic, data.tau);
  }
};

TDPCMSolver * createTDOnsagerIEFSolver(const TDSolverData & data) {
  buildTDOnsagerIEFSolver build;
  return for_id<derivative_types, integrator_types, TDPCMSolver>(
      build, data, data.howDerivative, data.howIntegrator);
}
const std::string TDONSAGERIEFSOLVER("TDONSAGERIEF");
const bool registeredTDOnsagerIEFSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDONSAGERIEFSOLVER,
        createTDOnsagerIEFSolver);
}

namespace {
struct buildTDSingleIEFSolver {
  template <typename T, typename U>
  TDPCMSolver * operator()(const TDSolverData & data) {
    return new TDSingleIEFSolver<T, U>(
        data.epsilonStatic, data.epsilonDynamic, data.tau, data.tauIEF);
  }
};

TDPCMSolver * createTDSingleIEFSolver(const TDSolverData & data) {
  buildTDSingleIEFSolver build;
  return for_id<derivative_types, integrator_types, TDPCMSolver>(
      build, data, data.howDerivative, data.howIntegrator);
}
const std::string TDSINGLEIEFSOLVER("TDSINGLEIEF");
const bool registeredTDSingleIEFSolver =
    Factory<TDPCMSolver, TDSolverData>::TheFactory().registerObject(
        TDSINGLEIEFSOLVER,
        createTDSingleIEFSolver);
}

#endif // REGISTERTDSOLVERTOFACTORY_HPP
