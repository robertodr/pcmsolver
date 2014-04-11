#ifndef METALSPHERE_HPP
#define METALSPHERE_HPP

#include <complex>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "IGreensFunction.hpp"

/*! \file MetalSphere.hpp
 *  \class MetalSphere
 *  \brief Class to describe spherical metal nanoparticles.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2011, 2014
 *
 *  This class is a wrapper around the Fortran routines written
 *  by Stefano Corni et al. to take into account the presence
 *  of metal nanoparticles.
 *  References:
 *  http://dx.doi.org/10.1063/1.1342241
 *  http://dx.doi.org/10.1063/1.1507579
 *  http://dx.doi.org/10.1063/1.1558036
 */

class MetalSphere : public GreensFunction<double>
{
private:
    typedef std::complex<double> dcomplex;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    MetalSphere(double eps, double epsRe, double epsIm,
                const Eigen::Vector3d & pos, double radius)
        : GreensFunction<double>(false), epsSolvent_(eps), epsMetal_(dcomplex(epsRe, epsIm)),
          sphPosition_(pos), sphRadius_(radius) {}
    virtual ~MetalSphere() {}
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual double compDiagonalElementS(double area) const ;
    virtual double compDiagonalElementD(double area, double radius) const;
    virtual double epsilon() const { return epsSolvent_; } // This is just to get it to compile...
    /*!
     *  Get the S and D matrices
     */
    virtual void operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const;
    /*!
     *  Get the S matrix
     */
    virtual void operator()(Eigen::MatrixXd & S,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas) const;
    /*!
     *  Get the D matrix
     */
    virtual void operator()(Eigen::MatrixXd & D,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const;

    friend std::ostream & operator<<(std::ostream & os, MetalSphere & gf) {
        return gf.printObject(os);
    }
private:
    virtual double evaluate(double * source, double * probe) const;
    double epsSolvent_;
    dcomplex epsMetal_;
    Eigen::Vector3d sphPosition_;
    double sphRadius_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    // The build functor and use of for_id are not necessary as MetalSphere
    // inherits from a GreensFunction<double>
    IGreensFunction * createMetalSphere(const greenData & _data)
    {
	// We pass some bogus arguments...
	Eigen::Vector3d orig;
	orig << 0.0, 0.0, 0.0;
        return new MetalSphere(_data.epsilon, 0.0, 0.0, orig, 1.0);
    }
    const std::string METALSPHERE("MetalSphere");
    const bool registeredMetalSphere =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            METALSPHERE, createMetalSphere);
}
#endif // METALSPHERE_HPP
