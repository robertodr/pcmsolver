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

#ifndef VCPCMSOLVER_HPP
#define VCPCMSOLVER_HPP

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;
class IGreensFunction;

#include "VPCMSolver.hpp"

/*! \file VCPCMSolver.hpp
 *  \class VCPCMSolver
 *  \brief A CPCM solver exploiting the variational formalism
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  We use the Non-Virtual Interface idiom.
 */

class VCPCMSolver __final : public VPCMSolver
{
public:
    VCPCMSolver(double corr) : VPCMSolver(), correction_(corr) {}
    virtual ~VCPCMSolver() {}

    friend std::ostream & operator<<(std::ostream & os, VCPCMSolver & solver) {
        return solver.printSolver(os);
    }
private:
    /*! Correction for the conductor results */
    double correction_;
    /*! CPCM S matrix, not symmetry blocked */
    Eigen::MatrixXd S_;
    /*! PCM matrix, symmetry blocked form */
    std::vector<Eigen::MatrixXd> blockS_;

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
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential, int irrep = 0) const __override attribute(const);
    virtual std::ostream & printSolver(std::ostream & os) __override;
};

#endif // VCPCMSOLVER_HPP
