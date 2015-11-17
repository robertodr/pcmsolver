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

#ifndef VIEFSOLVER_HPP
#define VIEFSOLVER_HPP

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class IGreensFunction;

#include "VPCMSolver.hpp"

/*! \file VIEFSolver.hpp
 *  \class VIEFSolver
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

class VIEFSolver __final : public VPCMSolver
{
public:
    VIEFSolver() : VPCMSolver() {}
    VIEFSolver(GuessType guess, UpdateType update) : VPCMSolver(guess, update) {}
    virtual ~VIEFSolver() {}

    /*! \brief Builds PCM matrix for an anisotropic environment
     *  \param[in] cavity the cavity to be used.
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    void buildAnisotropicMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o);
    /*! \brief Builds PCM matrix for an isotropic environment
     *  \param[in] cavity the cavity to be used.
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    void buildIsotropicMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o);

    friend std::ostream & operator<<(std::ostream & os, VIEFSolver & solver) {
        return solver.printSolver(os);
    }
private:
    /*! The \f$ \tilde{\mathbf{Y}} \f$ matrix */
    Eigen::MatrixXd tilde_Y_;
    /*! The \f$ \tilde{\mathbf{Y}} \f$ matrix in symmetry blocked form*/
    std::vector<Eigen::MatrixXd> blocktilde_Y_;
    /*! The \f$ \tilde{\mathbf{R}}_\infty \f$ matrix */
    Eigen::MatrixXd R_infinity_;
    /*! The \f$ \tilde{\mathbf{R}}_\infty \f$ matrix in symmetry blocked form */
    std::vector<Eigen::MatrixXd> blockR_infinity_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) __override;
    /*! \brief Updates the R^\dagger transformed ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd updateCharge_impl(const Eigen::VectorXd & potential, int irrep = 0) const __override attribute(const);
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] CGtol conjugate gradient solver tolerance
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential,
        double CGtol = Eigen::NumTraits<double>::epsilon(), int irrep = 0) const __override attribute(const);
    virtual std::ostream & printSolver(std::ostream & os) __override;

    /*! \brief A uniform ASC initial guess
     *  \param[in] nuc_chg total nuclear charge
     *  \param[in] irrep the irreducible representation of the ASC
     *  \return the initial guess vector
     *
     *  We suppose an initial uniform surface charge
     *  summing up to the total nuclear charge
     */
    virtual Eigen::VectorXd initialGuessUniform(double nuc_chg, int irrep = 0) const __override attribute(const);
    /*! \brief A diagonally scaled initial guess
     *  \param[in] potential the electrostatic potential
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     *  \return the initial guess vector
     *
     *  The initial guess for the ASC is calculated assuming a diagonal
     *  approximation for the PCM stiffness matrix
     */
    virtual Eigen::VectorXd initialGuessDiagonal(const Eigen::VectorXd & potential, int irrep = 0)
      const __override attribute(const);
    /*! \brief A low-accuracy initial guess
     *  \param[in] potential the electrostatic potential
     *  \param[in] irrep the irreducible representation
     *  \return the initial guess vector
     *
     *  The initial guess for the ASC is calculated with a low accuracy
     *  (10^-4) CG solver
     */
    virtual Eigen::VectorXd initialGuessLowAccuracy(const Eigen::VectorXd & potential, int irrep = 0)
      const __override attribute(const);

    /*! Transform ASC from the dressed to the bare representation */
    Eigen::VectorXd bareASC(const Eigen::VectorXd & dressedASC, int irrep = 0) const attribute(const) {
      int irrDim = blockR_infinity_[irrep].rows();
      Eigen::VectorXd bareASC = Eigen::VectorXd::Zero(dressedASC.size());
      bareASC.segment(irrep*irrDim, irrDim) = blockR_infinity_[irrep].adjoint() * dressedASC.segment(irrep*irrDim, irrDim);
      return bareASC;
    }
    /*! Transform MEP from the bare to the dressed representation */
    Eigen::VectorXd dressedMEP(const Eigen::VectorXd & bareMEP, int irrep = 0) const attribute(const) {
      int irrDim = blockR_infinity_[irrep].rows();
      Eigen::VectorXd dressedMEP = Eigen::VectorXd::Zero(bareMEP.size());
      dressedMEP.segment(irrep*irrDim, irrDim) = blockR_infinity_[irrep] * bareMEP.segment(irrep*irrDim, irrDim);
      return dressedMEP;
    }
};

#endif // VIEFSOLVER_HPP
