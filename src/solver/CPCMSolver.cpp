#include "CPCMSolver.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "Cavity.hpp"
#include "GreensFunction.hpp"

void CPCMSolver::buildSystemMatrix(Cavity & cavity) 
{
	if (greenInside->isUniform() && greenOutside->isUniform()) 
	{
		buildIsotropicMatrix(cavity);
	} 
	else 
	{
	        throw std::runtime_error("C-PCM is defined only for isotropic environments!");
	}
}


void CPCMSolver::buildIsotropicMatrix(Cavity & cav)
{
	double epsilon = greenOutside->getDielectricConstant();
    	int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	
	// This is the very core of PCMSolver
    	greenInside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
    	
	double fact = (epsilon - 1.0)/(epsilon + correction);
    	PCMMatrix = SI;
    	PCMMatrix = fact * PCMMatrix.inverse();
    	Eigen::MatrixXd PCMAdjoint(cavitySize, cavitySize); 
    	PCMAdjoint = PCMMatrix.adjoint().eval(); // See Eigen doc for the reason of this
    	PCMMatrix = 0.5 * (PCMMatrix + PCMAdjoint);
	// PRINT TO FILE RELEVANT INFO ABOUT PCMMatrix
    	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(PCMMatrix);
    	if (solver.info() != Eigen::Success) 
		abort();
    	std::ofstream matrixOut("PCM_matrix");
 	matrixOut << "PCM matrix printout" << std::endl;
    	matrixOut << "Number of Tesserae: " << cavitySize << std::endl;
    	matrixOut << "Largest Eigenvalue: " << solver.eigenvalues()[cavitySize-1] << std::endl;
    	matrixOut << "Lowest Eigenvalue: " << solver.eigenvalues()[0] << std::endl;
    	matrixOut << "Average of Eigenvalues: " << (solver.eigenvalues().sum() / cavitySize)<< std::endl;
    	matrixOut << "List of Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    	matrixOut.close();
    	builtIsotropicMatrix = true;
    	builtAnisotropicMatrix = false;
}

void CPCMSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) 
{
	if (builtIsotropicMatrix) 
	{
		charge = - PCMMatrix * potential;
	} 
	else 
	{
		throw std::runtime_error("PCM matrix not initialized!");
	}
}
    
std::ostream & CPCMSolver::printSolver(std::ostream & os) 
{
	os << "Solver Type: C-PCM";
	return os;
}

