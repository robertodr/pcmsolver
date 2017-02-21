#ifndef TDSOLVERHELPERFUNCTIONS_HPP
#define TDSOLVERHELPERFUNCTIONS_HPP

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

namespace td_solver {
inline Eigen::VectorXd K(const Eigen::VectorXd & Lambda, double factor) {
  Eigen::VectorXd K = Eigen::VectorXd::Zero(Lambda.size());
  for (int i = 0; i < Lambda.size(); ++i) {
    K(i) = (2 * M_PI - Lambda(i)) / (2 * M_PI * factor - Lambda(i));
  }
  return K;
}

inline Eigen::VectorXd tau(const Eigen::VectorXd & Lambda,
                           double e_d,
                           double e_0,
                           double tau_D) {
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(Lambda.size());
  for (int i = 0; i < Lambda.size(); ++i) {
    double num = (2 * M_PI - Lambda(i)) * e_d + 2 * M_PI + Lambda(i);
    double denom = (2 * M_PI - Lambda(i)) * e_0 + 2 * M_PI + Lambda(i);
    tau(i) = tau_D * num / denom;
  }
  return tau;
}

inline Eigen::VectorXd tauInverse(const Eigen::VectorXd & Lambda,
                                  double e_d,
                                  double e_0,
                                  double tau_D) {
  Eigen::VectorXd tau_inv = Eigen::VectorXd::Zero(Lambda.size());
  for (int i = 0; i < Lambda.size(); ++i) {
    double num = (2 * M_PI - Lambda(i)) * e_d + 2 * M_PI + Lambda(i);
    double denom = (2 * M_PI - Lambda(i)) * e_0 + 2 * M_PI + Lambda(i);
    tau_inv(i) = denom / (tau_D * num);
  }
  return tau_inv;
}
} // namespace td_solver

#endif // TDSOLVERHELPERFUNCTIONS_HPP
