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

    /*! Transform ASC from the dressed to the bare representation */
    Eigen::VectorXd getBareASC(const Eigen::VectorXd & dressedASC, int irrep = 0) const attribute(const);
    /*! Transform MEP from the bare to the dressed representation */
    Eigen::VectorXd getDressedMEP(const Eigen::VectorXd & bareMEP, int irrep = 0) const attribute(const);

    friend std::ostream & operator<<(std::ostream & os, VIEFSolver & solver) {
        return solver.printSolver(os);
    }
private:
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
    /*! \brief Calculates the error
     *  \param[in] dressedASC vector containing the dressed ASC at cavity points
     *  \param[in] bareMEP the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     *  \return error at the current iteration
     *
     *  This function calculates:
     *  \f[
     *    \mathbf{e}_Q^{(i)} = \tilde{\mathbf{Y}}\tilde{\mathbf{q}}^{(i)} + \tilde{\mathbf{v}}^{(i)} = -mathbf{r}^{(i)}
     *  \f]
     *  \note The error vector is the quantity to be used in the DIIS procedure.
     *  It is the residual vector but with opposite sign.
     */
    virtual Eigen::VectorXd error_impl(const Eigen::VectorXd & dressedASC,
        const Eigen::VectorXd & bareMEP, int irrep = 0) const __override;
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     *  \param[in] CGtol conjugate gradient solver tolerance
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential, int irrep = 0,
        double CGtol = Eigen::NumTraits<double>::epsilon()) const __override attribute(const);
    virtual std::ostream & printSolver(std::ostream & os) __override;

    /*! \brief Returns the initial guess for the ASC given the MEP and the desired irreducible representation
     *  \param[in] MEP the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd initialGuess_impl(const Eigen::VectorXd & MEP,
                                              double nuc_chg = 0.0, int irrep = 0) const __override;

    /*! \brief Updates the bare ASC given the dressed ASC, the residual and the desired irreducible representation
     *  \param[in] dressedASC vector containing the dressed ASC at cavity points
     *  \param[in] residual the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     *  \return the updated **bare** ASC
     *
     *  This function calculates the update
     *  \f[
     *    \tilde{\mathbf{q}}^{(i+1)} = \tilde{\mathbf{q}}^{(i)} + \alpha^{(i)}\mathbf{r}^{(i)}
     *  \f]
     *  and then transforms the updated ASC to the bare representation.
     *  \note This function takes care of the transformation to the bare representation
     *  of the ASC
     */
    virtual Eigen::VectorXd updateCharge_impl(const Eigen::VectorXd & dressedASC,
        const Eigen::VectorXd & residual, int irrep = 0) const __override;
};

#endif // VIEFSOLVER_HPP
