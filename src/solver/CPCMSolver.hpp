#ifndef CPCMSOLVER_HPP
#define CPCMSOLVER_HPP

#include <iosfwd>
#include <string>

#include <Eigen/Dense>

#include "Config.hpp"

class Cavity;
class GePolCavity;
class GreensFunction;

#include "PCMSolver.hpp"
#include "SolverFactory.hpp"

/*! \file CPCMSolver.hpp  
 *  \class CPCMSolver
 *  \brief Solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2013
 */

class CPCMSolver : public PCMSolver 
{
	private:
    		bool builtIsotropicMatrix;
    		bool builtAnisotropicMatrix;
    		double correction;
    		Eigen::MatrixXd PCMMatrix;
//    		static const double factor = 1.0694;
    		static const double factor = 1.07;
                void buildIsotropicMatrix(GePolCavity & cav);
    		virtual std::ostream & printSolver(std::ostream & os);
	public:
		CPCMSolver() {}
                CPCMSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, double correction_ = 0.0) 
			: PCMSolver(gfInside_, gfOutside_), builtIsotropicMatrix(false), builtAnisotropicMatrix(false), correction(correction_) {}                
                //CPCMSolver(const Section & solver);
                virtual ~CPCMSolver() {}
                const Eigen::MatrixXd & getPCMMatrix() const { return PCMMatrix; }
                virtual void buildSystemMatrix(Cavity & cavity);
                //virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                void setCorrection(double correction_) { correction = correction_; }
                friend std::ostream & operator<<(std::ostream & os, CPCMSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
};

namespace
{
	PCMSolver * createCPCMSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, double correction_ = 0.0, int integralEquation_ = 1)
	{
		return new CPCMSolver(gfInside_, gfOutside_, correction_);
	}
	const std::string CPCMSOLVER("CPCM");
	const bool registeredCPCMSolver = SolverFactory::TheSolverFactory().registerSolver(CPCMSOLVER, createCPCMSolver);
}

#endif // CPCMSOLVER_HPP
