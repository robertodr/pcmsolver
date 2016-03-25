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

#ifndef SOLVERIMPL_HPP
#define SOLVERIMPL_HPP

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"

/*! \file SolverImpl.cpp
 *  \brief Functions common to all solvers
 *  \author Roberto Di Remigio
 *  \date 2015
 */

/*! \brief Types of environment */
enum EnvironmentType {
  RegularIsotropic, /**< Regular isotropic environment: vacuum inside the cavity, uniform dielectric outside */
  FlippedIsotropic, /**< Flipped isotropic environment: vacuum outside the cavity, uniform dielectric inside */
  AnisotropicEnv    /**< Anisotropic environment: covers all other cases*/
};

inline EnvironmentType environment(const IGreensFunction & gf_i, const IGreensFunction & gf_o) {
  EnvironmentType retval = AnisotropicEnv;
  bool isotropic = (gf_i.uniform() && gf_o.uniform());
  // Determine if we are in the regular or flipped regime
  if (isotropic) {
    bool flipped = ((profiles::epsilon(gf_i.permittivity()) != 1.0) && (profiles::epsilon(gf_o.permittivity()) == 1.0));
    flipped ? retval = FlippedIsotropic : retval = RegularIsotropic;
  }
  return retval;
}

/*! \brief Builds the **anisotropic** IEFPCM matrix
 *  \param[in] cav the discretized cavity
 *  \param[in] gf_i Green's function inside the cavity
 *  \param[in] gf_o Green's function outside the cavity
 *  \return the \f$ \mathbf{K} = \mathbf{T}^{-1}\mathbf{R} \f$ matrix
 *
 *  This function calculates the PCM matrix. We use the following definitions:
 *  \f[
 *     \begin{align}
 *       \mathbf{T} &=
 *      \left(2\pi\mathbf{I} - \mathbf{D}_\mathrm{e}\mathbf{A}\right)\mathbf{S}_\mathrm{i}
 *      +\mathbf{S}_\mathrm{e}\left(2\pi\mathbf{I} + \mathbf{A}\mathbf{D}_\mathrm{i}^\dagger\right) \\
 *      \mathbf{R} &=
 *      \left(2\pi\mathbf{I} - \mathbf{D}_\mathrm{e} * \mathbf{A} \right) -
 *      \mathbf{S}_\mathrm{e}\mathbf{S}^{-1}_\mathrm{i}\left(2\pi\mathbf{I} - \mathbf{D}_\mathrm{i} * \mathbf{A}\right)
 *     \end{align}
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd anisotropicIEFMatrix(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
  TIMER_ON("Computing SI");
  Eigen::MatrixXd SI = gf_i.singleLayer(cav.elements());
  TIMER_OFF("Computing SI");
  TIMER_ON("Computing DI");
  Eigen::MatrixXd DI = gf_i.doubleLayer(cav.elements());
  TIMER_OFF("Computing DI");
  TIMER_ON("Computing SE");
  Eigen::MatrixXd SE = gf_o.singleLayer(cav.elements());
  TIMER_OFF("Computing SE");
  TIMER_ON("Computing DE");
  Eigen::MatrixXd DE = gf_o.doubleLayer(cav.elements());
  TIMER_OFF("Computing DE");

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(DE, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(SE, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }

  Eigen::MatrixXd A = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  // 1. Form T
  TIMER_ON("Assemble T matrix");
  Eigen::MatrixXd fullPCMMatrix = ((2 * M_PI * Id - DE * A) * SI + SE * (2 * M_PI * Id + A * DI.adjoint().eval()));
  TIMER_OFF("Assemble T matrix");
  // 2. Invert T using LU decomposition with full pivoting
  //    This is a rank-revealing LU decomposition, this allows us
  //    to test if T is invertible before attempting to invert it.
  TIMER_ON("Invert T matrix");
  Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix);
  if (!(T_LU.isInvertible())) PCMSOLVER_ERROR("T matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  fullPCMMatrix = T_LU.inverse();
  TIMER_OFF("Invert T matrix");
  Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
  if (!(SI_LU.isInvertible())) PCMSOLVER_ERROR("SI matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  TIMER_ON("Assemble T^-1R matrix");
  fullPCMMatrix *= ((2 * M_PI * Id - DE * A) - SE * SI_LU.inverse() * (2 * M_PI * Id - DI * A));
  TIMER_OFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}

/*! \brief Builds the **isotropic** IEFPCM matrix in a regular environment
 *  \param[in] cav the discretized cavity
 *  \param[in] gf  vacuum Green's function
 *  \param[in] epsilon permittivity outside the cavity
 *  \return the \f$ \mathbf{K} = \mathbf{T}^{-1}\mathbf{R} \f$ matrix
 *
 *  This function calculates the PCM matrix. We use the following definitions:
 *  \f[
 *     \begin{align}
 *       \mathbf{T} &=
 *      \left(2\pi\frac{\varepsilon+1}{\varepsilon-1}\mathbf{I} - \mathbf{D}\mathbf{A}\right)\mathbf{S} \\
 *      \mathbf{R} &= \left(2\pi\mathbf{I} - \mathbf{D} *\mathbf{A}\right)
 *     \end{align}
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd regularIsotropicIEFMatrix(const Cavity & cav, const IGreensFunction & gf, double epsilon)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute S and D on the whole cavity, regardless of symmetry
  TIMER_ON("Computing S");
  Eigen::MatrixXd S = gf.singleLayer(cav.elements());
  TIMER_OFF("Computing S");
  TIMER_ON("Computing D");
  Eigen::MatrixXd D = gf.doubleLayer(cav.elements());
  TIMER_OFF("Computing D");

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(D, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(S, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }

  Eigen::MatrixXd A = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  // Tq = -Rv -> q = -(T^-1 * R)v = -Kv
  // T = (2 * M_PI * fact * Id - D * A) * S; R = (2 * M_PI * Id - D * A)
  // fullPCMMatrix_ = K = T^-1 * R
  // 1. Form T
  double fact = (epsilon + 1.0)/(epsilon - 1.0);
  TIMER_ON("Assemble T matrix");
  Eigen::MatrixXd fullPCMMatrix = (2 * M_PI * fact * Id - D * A) * S;
  TIMER_OFF("Assemble T matrix");
  // 2. Invert T using LU decomposition with full pivoting
  //    This is a rank-revealing LU decomposition, this allows us
  //    to test if T is invertible before attempting to invert it.
  TIMER_ON("Invert T matrix");
  Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix);
  if (!(T_LU.isInvertible()))
    PCMSOLVER_ERROR("T matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  fullPCMMatrix = T_LU.inverse();
  TIMER_OFF("Invert T matrix");
  // 3. Multiply T^-1 and R
  TIMER_ON("Assemble T^-1R matrix");
  fullPCMMatrix *= (2 * M_PI * Id - D * A);
  TIMER_OFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}

/*! \brief Builds the **isotropic** IEFPCM matrix for a **flipped** environment
 *  \param[in] cav the discretized cavity
 *  \param[in] gf  vacuum Green's function
 *  \param[in] epsilon permittivity inside the cavity
 *  \return the \f$ \mathbf{K} = \mathbf{T}^{-1}\mathbf{R} \f$ matrix
 *
 *  This function calculates the PCM matrix. We use the following definitions:
 *  \f[
 *     \begin{align}
 *       \mathbf{T} &=
 *      \left(2\pi\frac{\varepsilon+1}{\varepsilon-1}\mathbf{I} + \mathbf{D}\mathbf{A}\right)\mathbf{S} \\
 *      \mathbf{R} &=
 *      -\left(2\pi\mathbf{I} - \mathbf{D}\mathbf{A}\right)
 *     \end{align}
 *  \f]
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd flippedIsotropicIEFMatrix(const Cavity & cav, const IGreensFunction & gf, double epsilon)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute S and D on the whole cavity, regardless of symmetry
  TIMER_ON("Computing S");
  Eigen::MatrixXd S = gf.singleLayer(cav.elements());
  TIMER_OFF("Computing S");
  TIMER_ON("Computing D");
  Eigen::MatrixXd D = gf.doubleLayer(cav.elements());
  TIMER_OFF("Computing D");

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(D, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(S, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }

  Eigen::MatrixXd A = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  // Tq = -Rv -> q = -(T^-1 * R)v = -Kv
  // T = (2 * M_PI * fact * Id + D * A) * S; R = -(2 * M_PI * Id - D * A)
  // fullPCMMatrix_ = K = T^-1 * R
  // 1. Form T
  double fact = (epsilon + 1.0)/(epsilon - 1.0);
  TIMER_ON("Assemble T matrix");
  Eigen::MatrixXd fullPCMMatrix = (2 * M_PI * fact * Id + D * A) * S;
  TIMER_OFF("Assemble T matrix");
  // 2. Invert T using LU decomposition with full pivoting
  //    This is a rank-revealing LU decomposition, this allows us
  //    to test if T is invertible before attempting to invert it.
  TIMER_ON("Invert T matrix");
  Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix);
  if (!(T_LU.isInvertible()))
    PCMSOLVER_ERROR("T matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  fullPCMMatrix = T_LU.inverse();
  TIMER_OFF("Invert T matrix");
  // 3. Multiply T^-1 and R
  TIMER_ON("Assemble T^-1R matrix");
  fullPCMMatrix *= -(2 * M_PI * Id - D * A);
  TIMER_OFF("Assemble T^-1R matrix");

  return fullPCMMatrix;
}

/*! \brief Builds the CPCM matrix
 *  \param[in] cav the discretized cavity
 *  \param[in] gf_i Green's function inside the cavity
 *  \param[in] epsilon permittivity outside the cavity
 *  \param[in] correction CPCM correction factor
 *  \return the \f$ \mathbf{K} = f(\varepsilon)\mathbf{S}^{-1} \f$ matrix
 *
 *  This function calculates the PCM matrix.
 *  The matrix is not symmetrized and is not symmetry packed.
 */
inline Eigen::MatrixXd CPCMMatrix(const Cavity & cav, const IGreensFunction & gf_i, double epsilon, double correction)
{
  // The total size of the cavity
  size_t cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();

  // Compute S on the whole cavity, regardless of symmetry
  TIMER_ON("Computing S");
  Eigen::MatrixXd S = gf_i.singleLayer(cav.elements());
  TIMER_OFF("Computing S");

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(S, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }

  double fact = (epsilon - 1.0)/(epsilon + correction);
  // Invert S  using LU decomposition with full pivoting
  // This is a rank-revealing LU decomposition, this allows us
  // to test if S is invertible before attempting to invert it.
  Eigen::FullPivLU<Eigen::MatrixXd> S_LU(S);
  if (!(S_LU.isInvertible()))
    PCMSOLVER_ERROR("S matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  return fact * S_LU.inverse();
}

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
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  // Form T
  return ((2 * M_PI * Id - DE * a) * SI + SE * (2 * M_PI * Id + a * DI.adjoint().eval()));
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

  // Compute S and D on the whole cavity, regardless of symmetry
  Eigen::MatrixXd S = gf_i.singleLayer(cav.elements());
  Eigen::MatrixXd D = gf_i.doubleLayer(cav.elements());

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the matrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    symmetryBlocking(D, cavitySize, dimBlock, nrBlocks);
    symmetryBlocking(S, cavitySize, dimBlock, nrBlocks);
  }

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  double fact = (epsilon + 1.0)/(epsilon - 1.0);
  return (2 * M_PI * fact * Id - D * a) * S;
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
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  // Form T
  Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
  if (!(SI_LU.isInvertible())) PCMSOLVER_ERROR("SI matrix is not invertible!", BOOST_CURRENT_FUNCTION);
  return ((2 * M_PI * Id - DE * a) - SE * SI_LU.inverse() * (2 * M_PI * Id - DI * a));
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

  // Compute D on the whole cavity, regardless of symmetry
  Eigen::MatrixXd D = gf_i.doubleLayer(cav.elements());

  // Perform symmetry blocking
  // If the group is C1 avoid symmetry blocking, we will just pack the matrix
  // into "block diagonal" when all other manipulations are done.
  if (cav.pointGroup().nrGenerators() != 0) {
    symmetryBlocking(D, cavitySize, dimBlock, nrBlocks);
  }

  Eigen::MatrixXd a = cav.elementArea().asDiagonal();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(cavitySize, cavitySize);

  return (2 * M_PI * Id - D * a);
}

#endif // SOLVERIMPL_HPP
