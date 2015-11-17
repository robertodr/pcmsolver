/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */


#ifndef VPCMSOLVER_HPP
#define VPCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

class Cavity;
class IGreensFunction;

/*! \file VPCMSolver.hpp
 *  \class VPCMSolver
 *  \brief Abstract Base Class for variational solvers inheritance hierarchy.
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  We use the Non-Virtual Interface idiom.
 *  "Bare" ASC and MEP are those **not** premultiplied by
 *  R_\infinity^\dagger, R_\infinity, respectively
 *  "Dressed" ASC and MEP are those premultiplied by
 *  R_\infinity^\dagger, R_\infinity, respectively
 */

class VPCMSolver
{
public:
    /*! \brief Initial ASC guess types */
    enum GuessType {
      Trivial, /**< Zero ASC */
      Uniform,  /**< Uniform ASC summing up to the solute nuclear charge */
      Diagonal, /**< ASC calculated from a diagonal approximation of the PCM matrix */
      LowAccuracy /**< ASC calculated from a low accuracy solution */
    };
    VPCMSolver() : built_(false), isotropic_(true), guess_(Trivial) {}
    VPCMSolver(GuessType guess) : built_(false), isotropic_(true), guess_(guess) {}
    virtual ~VPCMSolver() {}

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    void buildSystemMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) {
        buildSystemMatrix_impl(cavity, gf_i, gf_o);
    }
    /*! \brief Updates the R^\dagger transformed ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    Eigen::VectorXd updateCharge(const Eigen::VectorXd & potential, int irrep = 0) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet");
        return updateCharge_impl(potential, irrep);
    }
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    Eigen::VectorXd computeCharge(const Eigen::VectorXd & potential, int irrep = 0) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet");
        return computeCharge_impl(potential, irrep);
    }
    /*! \brief Returns the initial guess for the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    Eigen::VectorXd initialGuess(const Eigen::VectorXd & potential, double nuc_chg = 0.0, int irrep = 0) const {
        if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet");
        switch(guess_) {
          case Trivial:     return Eigen::VectorXd::Zero(potential.size());
          case Uniform:     return initialGuessUniform(nuc_chg, irrep);
          case Diagonal:    return initialGuessDiagonal(potential, irrep);
          case LowAccuracy: return initialGuessLowAccuracy(potential, irrep);
        }
    }

    friend std::ostream & operator<<(std::ostream & os, VPCMSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    /*! Whether the system matrix has been built */
    bool built_;
    /*! Whether the solver is isotropic */
    bool isotropic_;
    /*! Type of initial ASC guess */
    GuessType guess_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) = 0;
    /*! \brief Returns initial guess for the ASC
     */
    /*! \brief Updates the R^\dagger transformed ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd updateCharge_impl(const Eigen::VectorXd & potential, int irrep = 0) const = 0;
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] CGtol conjugate gradient solver tolerance
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential,
        double CGtol = Eigen::NumTraits<double>::epsilon(), int irrep = 0) const = 0;
    virtual std::ostream & printSolver(std::ostream & os) = 0;

    /*! \brief A uniform ASC initial guess
     *  \param[in] nuc_chg total nuclear charge
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     *  \return the initial guess vector
     *
     *  We suppose an initial uniform surface charge
     *  summing up to the total nuclear charge
     */
    virtual Eigen::VectorXd initialGuessUniform(double nuc_chg, int irrep = 0) const attribute(const) = 0;
    /*! \brief A diagonally scaled initial guess
     *  \param[in] potential the electrostatic potential
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     *  \return the initial guess vector
     *
     *  The initial guess for the ASC is calculated assuming a diagonal
     *  approximation for the PCM stiffness matrix
     */
    virtual Eigen::VectorXd initialGuessDiagonal(const Eigen::VectorXd & potential, int irrep = 0) const attribute(const) = 0;
    /*! \brief A low-accuracy initial guess
     *  \param[in] potential the electrostatic potential
     *  \param[in] irrep the irreducible representation
     *  \return the initial guess vector
     *
     *  The initial guess for the ASC is calculated with a low accuracy
     *  (10^-4) CG solver
     */
    virtual Eigen::VectorXd initialGuessLowAccuracy(const Eigen::VectorXd & potential, int irrep = 0) const attribute(const) = 0;
};

#endif // VPCMSOLVER_HPP
