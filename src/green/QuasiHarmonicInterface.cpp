#include "QuasiHarmonicInterface.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double QuasiHarmonicInterface<T>::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template<typename T>
T QuasiHarmonicInterface<T>::evaluate(T * sp, T * pp) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template <typename T>
void QuasiHarmonicInterface<T>::operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template <typename T>
void QuasiHarmonicInterface<T>::operator()(Eigen::MatrixXd & S,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template <typename T>
void QuasiHarmonicInterface<T>::operator()(Eigen::MatrixXd & D,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template<typename T>
double QuasiHarmonicInterface<T>::compDiagonalElementS(double area) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template<typename T>
double QuasiHarmonicInterface<T>::compDiagonalElementD(double area, double radius) const
{
    throw std::runtime_error("Green's function for a quasi-harmonic interface has not yet been implemented!");
}

template <typename T>
std::ostream & QuasiHarmonicInterface<T>::printObject(std::ostream & os)
{
    os << "Green's function type: quasi-harmonic interface" << std::endl;
    os << "Permittivity (layer 1) = " << eps1_ << std::endl;
    os << "Permittivity (layer 2) = " << eps2_ << std::endl;
    os << "Position               = " << pos_ << std::endl;
    os << "Width                  = " << width_;
    return os;
}

template class QuasiHarmonicInterface<double>;
template class QuasiHarmonicInterface<AD_directional>;
template class QuasiHarmonicInterface<AD_gradient>;
template class QuasiHarmonicInterface<AD_hessian>;
