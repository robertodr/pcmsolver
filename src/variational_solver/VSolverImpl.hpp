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

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "Cavity.hpp"
#include "Element.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"

/*! \file VSolverImpl.cpp
 *  \brief Functions common to all variational solvers
 *  \author Roberto Di Remigio
 *  \date 2015
 */

/*! \brief Builds the **anisotropic** \f$ \mathbf{T}_\varepsilon \f$ matrix
 *  \param[in] cav the discretized cavity
 *  \param[in] gf_i Green's function inside the cavity
 *  \param[in] gf_o Green's function outside the cavity
 *  \return the \f$ \mathbf{T}_\varepsilon \f$ matrix
 *
 *  We use the following definition:
 *  \f[
 *      \mathbf{T}_\varepsilon =
 *      \left(2\pi\mathbf{I} - \mathbf{D}_\mathrm{e}\mathbf{A}\right)\mathbf{S}_\mathrm{i}
 *      +\mathbf{S}_\mathrm{e}\left(2\pi\mathbf{I} +
 *      \mathbf{A}\mathbf{D}_\mathrm{i}^\dagger\right)
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd anisotropicTEpsilon(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  Eigen::MatrixXd SI = gf_i.singleLayer(cav.elements());
  Eigen::MatrixXd DI = gf_i.doubleLayer(cav.elements());
  Eigen::MatrixXd SE = gf_o.singleLayer(cav.elements());
  Eigen::MatrixXd DE = gf_o.doubleLayer(cav.elements());

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the matrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(DE, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SE, cavitySize, dimBlock, nrBlocks);
  }

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd aInv = a.inverse();

  // Form T
  return ((2 * M_PI * aInv - DE) * a * SI + SE * a * (2 * M_PI * aInv + DI.adjoint().eval()));
}

/*! \brief Builds the **isotropic** \f$ \mathbf{T}_\varepsilon \f$ matrix
 *  \param[in] cav the discretized cavity
 *  \param[in] gf_i Green's function inside the cavity
 *  \param[in] epsilon permittivity outside the cavity
 *  \return the \f$ \mathbf{T}_\varepsilon \f$ matrix
 *
 *  We use the following definition:
 *  \f[
 *      \mathbf{T}_\varepsilon =
 *      \left(2\pi\frac{\varepsilon+1}{\varepsilon-1}\mathbf{I} - \mathbf{D}_\mathrm{i}\mathbf{A}\right)\mathbf{S}_\mathrm{i}
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd isotropicTEpsilon(const Cavity & cav, const IGreensFunction & gf_i, double epsilon)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  Eigen::MatrixXd SI = gf_i.singleLayer(cav.elements());
  Eigen::MatrixXd DI = gf_i.doubleLayer(cav.elements());

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the matrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
  }

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd aInv = a.inverse();

  double fact = (epsilon + 1.0)/(epsilon - 1.0);
  return (2 * M_PI * fact * aInv - DI) * a * SI;
}

/*! \brief Builds the **anisotropic** \f$ \mathbf{R}_\infty \f$ matrix
 *  \param[in] cav the discretized cavity
 *  \param[in] gf_i Green's function inside the cavity
 *  \param[in] gf_o Green's function outside the cavity
 *  \return the \f$ \mathbf{R}_\infty\mathbf{A} \f$ matrix
 *
 *  We use the following definition:
 *  \f[
 *      \mathbf{R}_\infty =
 *      \left(2\pi\mathbf{A}^{-1} - \mathbf{D}_\mathrm{e}\right) -
 *      \mathbf{S}_\mathrm{e}\mathbf{S}^{-1}_\mathrm{i}\left(2\pi\mathbf{A}^{-1}-\mathbf{D}_\mathrm{i}\right)
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd anisotropicRinfinity(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  Eigen::MatrixXd SI = gf_i.singleLayer(cav.elements());
  Eigen::MatrixXd DI = gf_i.doubleLayer(cav.elements());
  Eigen::MatrixXd SE = gf_o.singleLayer(cav.elements());
  Eigen::MatrixXd DE = gf_o.doubleLayer(cav.elements());

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the matrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(DE, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SE, cavitySize, dimBlock, nrBlocks);
  }

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd aInv = a.inverse();

  // Form T
  Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
  if (!(SI_LU.isInvertible())) PCMSOLVER_ERROR("SI matrix is not invertible!");
  return (((2 * M_PI * aInv - DE) - SE * SI_LU.inverse() * (2 * M_PI * aInv - DI)) * a);
}

/*! \brief Builds the **isotropic** \f$ \mathbf{R}_\infty \f$ matrix
 *  \param[in] cav the discretized cavity
 *  \param[in] gf_i Green's function inside the cavity
 *  \return the \f$ \mathbf{R}_\infty\mathbf{A} \f$ matrix
 *
 *  We use the following definition:
 *  \f[
 *      \mathbf{R}_\infty =
 *      \left(2\pi\mathbf{A}^{-1} - \mathbf{D}_\mathrm{i}\right)
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd isotropicRinfinity(const Cavity & cav, const IGreensFunction & gf_i)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  Eigen::MatrixXd DI = gf_i.doubleLayer(cav.elements());

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the matrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
  }

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd aInv = a.inverse();

  return ((2 * M_PI * aInv - DI) * a);
}
