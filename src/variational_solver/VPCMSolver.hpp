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
    VPCMSolver() : built_(false), isotropic_(true) {}
    virtual ~VPCMSolver() {}

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    void buildSystemMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) {
        buildSystemMatrix_impl(cavity, gf_i, gf_o);
    }
    /*! \brief Returns initial guess for the ASC */
    Eigen::VectorXd initialGuess() const {
      if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet");
      return initialGuess_impl();
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

    friend std::ostream & operator<<(std::ostream & os, VPCMSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    /*! Whether the system matrix has been built */
    bool built_;
    /*! Whether the solver is isotropic */
    bool isotropic_;

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
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential, int irrep = 0) const = 0;
    virtual std::ostream & printSolver(std::ostream & os) = 0;
};

#endif // VPCMSOLVER_HPP
