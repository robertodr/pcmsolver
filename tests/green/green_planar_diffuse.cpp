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

#include "catch.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>


#include <Eigen/Core>

#include "bi_operators/CollocationIntegrator.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/PlanarDiffuse.hpp"
#include "green/dielectric_profile/OneLayerErf.hpp"
#include "green/dielectric_profile/OneLayerTanh.hpp"

SCENARIO("Evaluation of the planar diffuse Green's function and its derivatives", "[green][green_planar_diffuse]")
{
    GIVEN("A permittivity profile modelled by the hyperbolic tangent function")
    {
        // High dielectric constant inside
        double eps1 = 80.0;
        // Low dielectric constant outside
        double eps2 = 2.0;
        double width = 15.0;
        // Evaluation at negative z
        Eigen::Vector3d source1 = (Eigen::Vector3d() << 1.0, 0.0, -3.0).finished();
        Eigen::Vector3d sourceNormal1 = source1;
        sourceNormal1.normalize();
        Eigen::Vector3d probe1 = (Eigen::Vector3d() << 2.0, 0.0, 3.0).finished();
        Eigen::Vector3d probeNormal1 = probe1;
        probeNormal1.normalize();
        // Evaluation at positive z
        Eigen::Vector3d source2 = (Eigen::Vector3d() << 1.0, 0.0, 3.0).finished();
        Eigen::Vector3d sourceNormal2 = source2;
        sourceNormal2.normalize();
        Eigen::Vector3d probe2 = (Eigen::Vector3d() << 2.0, 0.0, 3.0).finished();
        Eigen::Vector3d probeNormal2 = probe2;
        probeNormal2.normalize();
        WHEN("the plane is at z=0")
        {
            double zInt = 0.0;
            PlanarDiffuse<> gf(eps1, eps2, width, zInt);
            THEN("the value of the Green's function inside the droplet is")
            {
                double value = 0.012507311388168523;
                double gf_value = gf.kernelS(source1, probe1);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            AND_THEN("the value of the Green's function outside the droplet is")
            {
                double value = 0.50004567416494572;
                double gf_value = gf.kernelS(source2, probe2);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            THEN("the value of the Green's function directional derivative wrt the probe point inside the droplet is")
            {
                double derProbe = -0.012506305835640469;
                double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the probe point outside the droplet is")
            {
                double derProbe = -0.29005549287308696;
                double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            THEN("the value of the Green's function directional derivative wrt the source point inside the droplet is")
            {
                double derSource = 0.012498621932118328;
                double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the source point outside the droplet is")
            {
                double derSource = 0.28879776627577236;
                double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
        }
		/*
        AND_WHEN("the plane is not at z=0")
        {
			double zInt = 10.0;
            PlanarDiffuse<> gf(eps1, eps2, width, zInt);
            THEN("the value of the Green's function inside the droplet is")
            {
                double value = 0.0125233347669694017;
                double gf_value = gf.kernelS(source1, probe1);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            AND_THEN("the value of the Green's function outside the droplet is")
            {
                double value = 0.5000329900631173;
                double gf_value = gf.kernelS(source2, probe2);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            THEN("the value of the Green's function directional derivative wrt the probe point inside the droplet is")
            {
                double derProbe = -0.0125024363466751109;
                double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the probe point outside the droplet is")
            {
                double derProbe = -0.289999709125188243;
                double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            THEN("the value of the Green's function directional derivative wrt the source point inside the droplet is")
            {
                double derSource = 0.0124899030052444404;
                double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the source point outside the droplet is")
            {
                double derSource = 0.288657584958107449;
                double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
        }
    }

    GIVEN("A permittivity profile modelled by the error function")
    {
        // High dielectric constant inside
        double eps1 = 80.0;
        // Low dielectric constant outside
        double eps2 = 2.0;
        double width = 5.0;
        // Evaluation at negative z
        Eigen::Vector3d source1 = (Eigen::Vector3d() << 1.0, 0.0, -3.0).finished();
        Eigen::Vector3d sourceNormal1 = source1;
        sourceNormal1.normalize();
        Eigen::Vector3d probe1 = (Eigen::Vector3d() << 2.0, 0.0, 3.0).finished();
        Eigen::Vector3d probeNormal1 = probe1;
        probeNormal1.normalize();
        // Evaluation at positive z
        Eigen::Vector3d source2 = (Eigen::Vector3d() << 1.0, 0.0, 3.0).finished();
        Eigen::Vector3d sourceNormal2 = source2;
        sourceNormal2.normalize();
        Eigen::Vector3d probe2 = (Eigen::Vector3d() << 2.0, 0.0, 3.0).finished();
        Eigen::Vector3d probeNormal2 = probe2;
        probeNormal2.normalize();
        WHEN("the plane is at z=0")
        {
			double zInt = 0.0;
            PlanarDiffuse<CollocationIntegrator, OneLayerErf> gf(eps1, eps2, width, zInt);
            THEN("the value of the Green's function inside the droplet is")
            {
                double value = 0.012507311377769814;
                double gf_value = gf.kernelS(source1, probe1);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            AND_THEN("the value of the Green's function outside the droplet is")
            {
                double value = 0.49991229576650942;
                double gf_value = gf.kernelS(source2, probe2);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            THEN("the value of the Green's function directional derivative wrt the probe point inside the droplet is")
            {
                double derProbe = -0.012506340818360315;
                double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the probe point outside the droplet is")
            {
                double derProbe = -0.28997553534915177;
                double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            THEN("the value of the Green's function directional derivative wrt the source point inside the droplet is")
            {
                double derSource = 0.012498621813619368;
                double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the source point outside the droplet is")
            {
                double derSource = 0.28871292628573908;
                double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
        }

        AND_WHEN("the plane is not at z=0")
        {
			double zInt = 10.0;
				PlanarDiffuse<CollocationIntegrator, OneLayerErf> gf(eps1, eps2, width, zInt);
            THEN("the value of the Green's function inside the droplet is")
            {
                double value = 0.012523344896520634;
                double gf_value = gf.kernelS(source1, probe1);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            AND_THEN("the value of the Green's function outside the droplet is")
            {
                double value = 0.49989736527661349;
                double gf_value = gf.kernelS(source2, probe2);
                INFO("ref_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
                INFO("gf_value  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
                REQUIRE(value == Approx(gf_value));
            }
            THEN("the value of the Green's function directional derivative wrt the probe point inside the droplet is")
            {
                double derProbe = -0.012502439280500169;
                double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the probe point outside the droplet is")
            {
                double derProbe = -0.28995345687399254;
                double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
                INFO("ref_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
                INFO("gf_derProbe  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
                REQUIRE(derProbe == Approx(gf_derProbe));
            }
            THEN("the value of the Green's function directional derivative wrt the source point inside the droplet is")
            {
                double derSource = 0.012489903372303254;
                double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
            AND_THEN("the value of the Green's function directional derivative wrt the source point outside the droplet is")
            {
                double derSource = 0.28869402146192158;
                double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
                INFO("ref_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
                INFO("gf_derSource  = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
                REQUIRE(derSource == Approx(gf_derSource));
            }
        }
		*/
    }
}
