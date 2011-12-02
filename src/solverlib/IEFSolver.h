/*! \file IEFSolver.h 
\brief PCM solver
*/


#ifndef IEFSOLVER
#define IEFSOLVER

#include <string>
#include <vector>
#include <iostream>
#include <complex>

using namespace std;

class IEFSolver : public PCMSolver {
 public:
    IEFSolver(GreensFunction &gfi, GreensFunction &gfo);
    IEFSolver(GreensFunction *gfi, GreensFunction *gfo);
    IEFSolver(Section solver);
    ~IEFSolver();
    const MatrixXd& getPCMMatrix() const {return PCMMatrix;};

    virtual void buildSystemMatrix(Cavity & cavity);
    virtual void buildAnisotropicMatrix(GePolCavity cav);
    virtual void buildIsotropicMatrix(GePolCavity cav);
    virtual double compDiagonalElementSoper(GreensFunction *green, int i, 
                                            GePolCavity cav);
    virtual double compDiagonalElementDoper(GreensFunction *green, int i, 
                                            GePolCavity cav);
    virtual VectorXd compCharge(const VectorXd & potential);
    virtual void compCharge(const VectorXd & potential, VectorXd & charge);
 private:
    bool builtAnisotropicMatrix;
    bool builtIsotropicMatrix;
    //    static const double factor = 1.0694;
    static const double factor = 1.07;
    MatrixXd PCMMatrix;
};
#endif
