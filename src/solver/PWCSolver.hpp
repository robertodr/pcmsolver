/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef PWCSOLVER_HPP
#define PWCSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "SparseMatrix.hpp"
#include "GenericAnsatzFunction.hpp"
#include "ConAnsatzFunction.hpp"
#include "readPoints.hpp"

class IGreensFunction;
class Vector2;
class Vector3;
class Cavity;
class WaveletCavity;
struct Compression;

#include "PCMSolver.hpp"

/*! \file PWCSolver.hpp
 *  \class PWCSolver
 *  \brief Class describing a wavelet solver with piecewise constant wavelet functions
 *  \author Luca Frediani and Monica Bugeanu
 *  \date 2014
 */

class PWCSolver : public PCMSolver
{
private:
    unsigned int interpolationGrade;
    unsigned int interpolationType;
    void initWEMMembers();
    virtual std::ostream & printSolver(std::ostream & os) __override;
public:
    PWCSolver(const Compression & _comp, int integralEquation_ = SecondKind)
        : PCMSolver(), interpolationGrade(2), interpolationType(1),
	af( new ConAnsatzFunction(_comp) ), integralEquation(integralEquation_) {
        initWEMMembers();
    }
    PWCSolver(const PWCSolver& old);
    virtual ~PWCSolver();
    Interpolation* getT_() {
        return af->interCoeff;
    }
    int getQuadratureLevel() {
        return af->quadratureLevel_;
    }
    friend std::ostream & operator<<(std::ostream & os, PWCSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    virtual void buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) __override;
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential, int irrep = 0) const __override;
    void constructSystemMatrix(const IGreensFunction & gf_i, const IGreensFunction & gf_o);
    void uploadCavity(const WaveletCavity & cavity); // different interpolation
    void constructSi(const IGreensFunction & gf_i);
    void constructSe(const IGreensFunction & gf_o);
    void solveFirstKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
    void solveSecondKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
    void solveFull(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
    void constructWavelets();
    void initInterpolation();

    // This is a dirty hack!!
    double epsilon_;
    double threshold;
    SparseMatrix S_i_, S_e_; // System matrices
    bool systemMatricesInitialized_;
    GenericAnsatzFunction * af;

    Vector3 *** pointList; // the old U

    unsigned int apriori1_, aposteriori1_;    // System matrix sparsities
    unsigned int apriori2_, aposteriori2_;    // System matrix sparsities
    int integralEquation;

    enum EquationType {FirstKind, SecondKind, Full};
};

#endif // PWCSOLVER_HPP
