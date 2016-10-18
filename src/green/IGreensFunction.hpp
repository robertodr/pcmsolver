/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#ifndef IGREENSFUNCTION_HPP
#define IGREENSFUNCTION_HPP

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Element;

#include "dielectric_profile/ProfileTypes.hpp"

/*! \file IGreensFunction.hpp
 *  \class IGreensFunction
 *  \brief Interface for Green's function classes
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2016
 *
 *  We **define** as _Green's function_ a function:
 *  \f[
 *      G(\mathbf{r}, \mathbf{r}^\prime) : \mathbb{R}^6 \rightarrow \mathbb{R}
 *  \f]
 *  that can be written as the sum of a singular and a non-singular part:
 *  \f[
 *      G(\mathbf{r}, \mathbf{r}^\prime) =
 *      \frac{1}{C(\mathbf{r}, \mathbf{r}^\prime)|\mathbf{r} - \mathbf{r}^\prime|} +
 *      G_\mathrm{img}(\mathbf{r}, \mathbf{r}^\prime)
 *  \f]
 *  The singular part is the Coulomb portion of the Green's function, possibly
 *  scaled by a position-dependent scalar coefficient:
 *  \f$C(\mathbf{r}, \mathbf{r}^\prime)\f$
 *  The non-singular part is usually associated with the image portion of the
 *  Green's function. This does not necessarily mean that the non-singular
 *  portion is physically an image contribution.
 *
 *  # Design motivation and rationale
 *
 *  Our definition of an electrostatic Green's function is rather arbitrary.
 *  However, at the time of writing, it covers all possible electrostatic
 *  Green's function already available within the library.
 *  Why do we need such a design?
 *  Green's functions and their directional derivatives appear as kernels of
 *  the \f$\mathcal{S}\f$ and \f$\mathcal{D}\f$ integral operators.
 *  Forming the matrix representation of these operators requires performing an
 *  integration over a surface finite element.
 *  Since these Green's functions present a Coulombic divergence, the diagonal
 *  elements of the operators will diverge unless appropriately formulated.
 *  This is possible, but requires **explicit** access to the expressions for
 *  the Coulomb singularity scaling coefficient and non-singular part. Hiding
 *  these details would otherwise introduce tight coupling with the solver
 *  and/or an intermediate boundary integral operator object.
 *  The Non-Virtual Interface (NVI) idiom is used.
 */

/*! \typedef KernelS
 *  \brief functor handle to the kernelS method
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
    KernelS;

/*! \typedef KernelD
 *  \brief functor handle to the kernelD method
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &,
                             const Eigen::Vector3d &)> KernelD;

class IGreensFunction {
public:
  virtual ~IGreensFunction() {}

  /**@{ Methods to sample the Green's function and its probe point directional
   * derivative */
  /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e.
   * the value of the
   *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1,
   * \mathbf{p}_2)\f$
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double kernelS(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
    return kernelS_impl(p1, p2);
  }
  /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the
   * pair of points p1, p2:
   *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1,
   * \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
   *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this
   * methods with \f$\mathbf{p}_1\f$
   *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} =
   * \mathbf{n}_{\mathbf{p}_1}\f$
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double kernelD(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1,
                 const Eigen::Vector3d & p2) const {
    return kernelD_impl(direction, p1, p2);
  }
  /**@}*/

  KernelS exportKernelS() const { return exportKernelS_impl(); }
  KernelD exportKernelD() const { return exportKernelD_impl(); }

  /**@{ Methods to sample the Green's function singular part and its probe point
   * directional derivative */
  /*! \brief Returns Coulomb singularity separation coefficient
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double coefficientCoulomb(const Eigen::Vector3d & source,
                            const Eigen::Vector3d & probe) const {
    return coefficient_impl(source, probe);
  }
  /*! \brief Returns singular part of the Green's function
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   *  \note This method only has a Non-Virtual Interface (NVI)
   */
  double Coulomb(const Eigen::Vector3d & source,
                 const Eigen::Vector3d & probe) const {
    double r12 = (source - probe).norm();
    return (1.0 / (coefficient_impl(source, probe) * r12));
  }
  /*! Returns value of the directional derivative of the
   *  Coulomb singularity separation coefficient for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double coefficientCoulombDerivative(const Eigen::Vector3d & direction,
                                      const Eigen::Vector3d & p1,
                                      const Eigen::Vector3d & p2) const {
    return coefficientCoulombDerivative_impl(direction, p1, p2);
  }
  /*! Returns value of the directional derivative of the
   *  singular part of the Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double CoulombDerivative(const Eigen::Vector3d & direction,
                           const Eigen::Vector3d & p1,
                           const Eigen::Vector3d & p2) const {
    return CoulombDerivative_impl(direction, p1, p2);
  }
  /**@}*/

  /**@{ Methods to sample the Green's function non-singular part (aka the image) and
   * its probe point directional derivative */
  /*! \brief Returns non-singular part of the Green's function (image potential)
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double imagePotential(const Eigen::Vector3d & source,
                        const Eigen::Vector3d & probe) const {
    return imagePotential_impl(source, probe);
  }
  /*! Returns value of the directional derivative of the
   *  non-singular part (image potential) of the Greens's function for the pair of
   * points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double imagePotentialDerivative(const Eigen::Vector3d & direction,
                                  const Eigen::Vector3d & p1,
                                  const Eigen::Vector3d & p2) const {
    return imagePotentialDerivative_impl(direction, p1, p2);
  }
  /**@}*/

  /*! Whether the Green's function describes a uniform environment */
  virtual bool uniform() const = 0;
  /*! Returns a dielectric permittivity profile */
  virtual Permittivity permittivity() const = 0;

  friend std::ostream & operator<<(std::ostream & os, IGreensFunction & gf) {
    return gf.printObject(os);
  }

protected:
  /**@{ Methods to sample the Green's function and its probe point directional
   * derivative */
  /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e.
   * the value of the
   *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1,
   * \mathbf{p}_2)\f$
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   */
  virtual double kernelS_impl(const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const = 0;
  /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the
   * pair of points p1, p2:
   *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1,
   * \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
   *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this
   * methods with \f$\mathbf{p}_1\f$
   *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} =
   * \mathbf{n}_{\mathbf{p}_1}\f$
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const = 0;
  /**@}*/

  virtual KernelS exportKernelS_impl() const = 0;
  virtual KernelD exportKernelD_impl() const = 0;

  /**@{ Methods to sample the Green's function singular part and its probe point
   * directional derivative */
  /*! \brief Returns Coulomb singularity separation coefficient
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   */
  virtual double coefficient_impl(const Eigen::Vector3d & source,
                                  const Eigen::Vector3d & probe) const = 0;
  /*! Returns value of the directional derivative of the
   *  Coulomb singularity separation coefficient for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double coefficientCoulombDerivative_impl(
      const Eigen::Vector3d & direction, const Eigen::Vector3d & p1,
      const Eigen::Vector3d & p2) const = 0;
  /*! Returns value of the directional derivative of the
   *  singular part of the Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double CoulombDerivative_impl(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1,
                                        const Eigen::Vector3d & p2) const = 0;
  /**@}*/

  /**@{ Methods to sample the Green's function non-singular part (aka the image) and
   * its probe point directional derivative */
  /*! \brief Returns non-singular part of the Green's function (image potential)
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   */
  virtual double imagePotential_impl(const Eigen::Vector3d & source,
                                     const Eigen::Vector3d & probe) const = 0;
  /*! Returns value of the directional derivative of the
   *  non-singular part (image potential) of the Greens's function for the pair of
   * points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  virtual double imagePotentialDerivative_impl(const Eigen::Vector3d & direction,
                                               const Eigen::Vector3d & p1,
                                               const Eigen::Vector3d & p2) const = 0;
  /**@}*/
  virtual std::ostream & printObject(std::ostream & os) = 0;
};

#endif // IGREENSFUNCTION_HPP
