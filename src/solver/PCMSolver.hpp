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

#ifndef PCMSOLVER_HPP
#define PCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"

class Cavity;
class IGreensFunction;

#include "SolverImpl.hpp"

/*! \file PCMSolver.hpp
 *  \class PCMSolver
 *  \brief Abstract Base Class for solvers inheritance hierarchy.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2011, 2015
 *
 *  We use the Non-Virtual Interface idiom.
 */

class PCMSolver
{
  public:
    /*! \note Initialization takes the worst possible scenario into account:
     *  system matrix not build and anisotropic environment
     */
    PCMSolver() : built_(false), environment_(AnisotropicEnv) {}
    virtual ~PCMSolver() {}

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     *  \note We initialize environment_ to the correct value
     */
    void buildSystemMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) {
      // Determine which type of environment we will have to handle
      environment_ = environment(gf_i, gf_o);
      buildSystemMatrix_impl(cavity, gf_i, gf_o);
    }
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    Eigen::VectorXd computeCharge(const Eigen::VectorXd & potential, int irrep = 0) const {
      if (!built_) PCMSOLVER_ERROR("PCM matrix not calculated yet", BOOST_CURRENT_FUNCTION);
      return computeCharge_impl(potential, irrep);
    }

    friend std::ostream & operator<<(std::ostream & os, PCMSolver & solver) {
      return solver.printSolver(os);
    }
  protected:
    /*! Whether the system matrix has been built */
    bool built_;
    /*! Type of environment */
    EnvironmentType environment_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) = 0;
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential, int irrep = 0) const = 0;
    virtual std::ostream & printSolver(std::ostream & os) = 0;
};

#endif // PCMSOLVER_HPP
