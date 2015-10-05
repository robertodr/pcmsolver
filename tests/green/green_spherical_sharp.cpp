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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "AnalyticEvaluate.hpp"
#include "SphericalSharp.hpp"
#include "DerivativeTypes.hpp"

TEST_CASE("Evaluation of the spherically symmetric sharp interface Green's function and its derivatives", "[green][green_spherical_sharp]")
{
    double epsilonSolvent = 1.0;
    double eps = 78.39;
    double radius = 50.0;
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d source = 100 * Eigen::Vector3d::Random();
    Eigen::Vector3d sourceNormal = source + Eigen::Vector3d::Random();
    sourceNormal.normalize();
    Eigen::Vector3d probe = 100 * Eigen::Vector3d::Random();
    Eigen::Vector3d probeNormal = probe + Eigen::Vector3d::Random();
    probeNormal.normalize();
    Eigen::Array4d result = analyticSphericalSharp(eps, epsilonSolvent, radius, origin, sourceNormal, source, probeNormal, probe);

    /*! \class SphericalSharp
     *  \test \b SphericalSharpTest_numerical tests the numerical evaluation of the SphericalSharp Green's function against analytical result
     */
    SECTION("Numerical derivative")
    {
        SphericalSharp<Numerical, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
        double analytic = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(analytic == Approx(gf_value));

        double analytic_derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(analytic_derProbe == Approx(gf_derProbe));

        double analytic_derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(analytic_derSource == Approx(gf_derSource));
    }

    /*! \class SphericalSharp
     *  \test \b SphericalSharpTest_directional_AD tests the automatic evaluation (directional derivative only)
     *  of the SphericalSharp Green's function against analytical result
     */
    SECTION("Directional derivative via AD")
    {
        SphericalSharp<AD_directional, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
        double analytic = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(analytic == Approx(gf_value));

        double analytic_derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(analytic_derProbe == Approx(gf_derProbe));

        double analytic_derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(analytic_derSource == Approx(gf_derSource));
    }

    /*! \class SphericalSharp
     *  \test \b SphericalSharpTest_gradient_AD tests the automatic evaluation (full gradient)
     *  of the SphericalSharp Green's function against analytical result
     */
    SECTION("Gradient via AD")
    {
        SphericalSharp<AD_gradient, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
        double analytic = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(analytic == Approx(gf_value));

        double analytic_derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(analytic_derProbe == Approx(gf_derProbe));

        double analytic_derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(analytic_derSource == Approx(gf_derSource));
    }

    /*! \class SphericalSharp
     *  \test \b SphericalSharpTest_hessian_AD tests the automatic evaluation (full hessian)
     *  of the SphericalSharp Green's function against analytical result
     */
    SECTION("Hessian via AD")
    {
        SphericalSharp<AD_hessian, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
        double analytic = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(analytic == Approx(gf_value));

        double analytic_derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(analytic_derProbe == Approx(gf_derProbe));

        double analytic_derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(analytic_derSource == Approx(gf_derSource));
    }
}
