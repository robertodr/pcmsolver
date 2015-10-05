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
#include "MetalNP.hpp"
#include "DerivativeTypes.hpp"

TEST_CASE("Evaluation of the spherical metal nanoparticle Green's function and its derivatives", "[green][green_metalnp]")
{
    double epsilonSolvent = 80.0;
    double epsRe = -9.906;
    double epsIm = -4.329;
    double radius = 50.0;
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d source = Eigen::Vector3d::Random();
    Eigen::Vector3d sourceNormal = source + Eigen::Vector3d::Random();
    sourceNormal.normalize();
    Eigen::Vector3d probe = Eigen::Vector3d::Random();
    Eigen::Vector3d probeNormal = probe + Eigen::Vector3d::Random();
    probeNormal.normalize();
    Eigen::Array4d result = analyticUniformDielectric(epsilonSolvent, sourceNormal, source, probeNormal, probe);
    //Eigen::Array4d result = analyticMetalNP(epsilonSolvent, espRe, epsIm, sourceNormal, source, probeNormal, probe);

    /*! \class MetalNP
     *  \test \b MetalNPTest_numerical tests the numerical evaluation of the MetalNP Green's function against analytical result
     */
    SECTION("Numerical derivative")
    {
        MetalNP<Numerical, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));

        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));

        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
    }

    /*! \class MetalNP
     *  \test \b MetalNPTest_directional_AD tests the automatic evaluation (directional derivative only)
     *  of the MetalNP Green's function against analytical result
     */
    SECTION("Directional derivative via AD")
    {
        MetalNP<AD_directional, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));

        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));

        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
    }

    /*! \class MetalNP
     *  \test \b MetalNPTest_gradient_AD tests the automatic evaluation (full gradient)
     *  of the MetalNP Green's function against analytical result
     */
    SECTION("Gradient via AD")
    {
        MetalNP<AD_gradient, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));

        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));

        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
    }

    /*! \class MetalNP
     *  \test \b MetalNPTest_hessian_AD tests the automatic evaluation (full hessian)
     *  of the MetalNP Green's function against analytical result
     */
    SECTION("Hessian via AD")
    {
        MetalNP<AD_hessian, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));

        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));

        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
    }
}

