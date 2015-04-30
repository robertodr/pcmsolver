#define BOOST_TEST_MODULE GreensFunctionTanhSphericalDiffuse

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "TanhSphericalDiffuse.hpp"

struct TanhSphericalDiffuseTest {
    double epsInside, epsOutside, sphereRadius, width;
    double inside_reference, outside_reference;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    TanhSphericalDiffuseTest() { SetUp(); }
    void SetUp() {
        epsInside = 80.0;
        epsOutside = 2.0;
	    sphereCenter << 0.0, 0.0, 0.0;
        sphereRadius = 100.0;
	    width = 5.0;
	    // Evaluation inside the sphere
        source1 << 1.0, 0.0, 0.0;
        sourceNormal1 = source1; // + Eigen::Vector3d::Random();
        sourceNormal1.normalize();
        probe1 << 2.0, 0.0, 0.0;
        probeNormal1 = probe1; // + Eigen::Vector3d::Random();
        probeNormal1.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        inside_reference = 0.012377848015483936;
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2; // + Eigen::Vector3d::Random();
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2; // + Eigen::Vector3d::Random();
        probeNormal2.normalize();
        // Reference value
        // Checked by comparing the asymptotic behaviour
        outside_reference = 0.50008829802731714;
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuse1, TanhSphericalDiffuseTest)

/*! \class TanhSphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseTest)
{
    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);
    double value = inside_reference;
    double gf_value = gf.function(source1, probe1);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_CHECK_CLOSE(gf_value, value, 1.0e-08);

    /*
    double derProbe = resultInside(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("derProbe    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_CHECK_CLOSE(derProbe, gf_derProbe, 1.0e-06);

    double derSource = resultInside(2);
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("derSource    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_CHECK_CLOSE(derSource, gf_derSource, 1.0e-06);
    */
}
/*! \class TanhSphericalDiffuse
 *  \test \b TanhSphericalDiffuseTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseTest)
{
    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);
    double value = outside_reference;
    double gf_value = gf.function(source2, probe2);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_CHECK_CLOSE(gf_value, value, 1.0e-08);

    /*
    double derProbe = resultOutside(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("derProbe    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_CHECK_CLOSE(derProbe, gf_derProbe, 1.0e-06);

    double derSource = resultOutside(2);
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("derSource    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_CHECK_CLOSE(derSource, gf_derSource, 1.0e-06);
    */
}
BOOST_AUTO_TEST_SUITE_END()

struct TanhSphericalDiffuseBogusTest {
    Eigen::Array4d analyticEvaluate(double eps, const Eigen::Vector3d & spNormal,
                                    const Eigen::Vector3d & sp,
                                    const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
        Eigen::Array4d result = Eigen::Array4d::Zero();
        double distance = (sp - pp).norm();
        double distance_3 = std::pow(distance, 3);
        double distance_5 = std::pow(distance, 5);

        // Value of the function
        result(0) = 1.0 / (eps * distance);
        // Value of the directional derivative wrt probe
        result(1) = (sp - pp).dot(ppNormal) / (eps * distance_3);
        // Directional derivative wrt source
        result(2) = - (sp - pp).dot(spNormal) / (eps * distance_3);
        // Value of the Hessian
        result(3) = spNormal.dot(ppNormal) / (eps * distance_3) - 3 * ((
                        sp - pp).dot(spNormal))*((sp - pp).dot(
                                    ppNormal)) / (eps * distance_5);

        return result;
    }
    double epsInside, epsOutside, sphereRadius, width;
    Eigen::Vector3d sphereCenter;
    Eigen::Vector3d source1, probe1, sourceNormal1, probeNormal1;
    Eigen::Vector3d source2, probe2, sourceNormal2, probeNormal2;
    Eigen::Array4d resultInside, resultOutside;
    TanhSphericalDiffuseBogusTest() { SetUp(); }
    void SetUp() {
        epsInside = 80.0;
        epsOutside = 80.0;
	    sphereCenter << 0.0, 0.0, 0.0;
        sphereRadius = 100.0;
	    width = 5.0;
	    // Evaluation inside the sphere
        source1 << 1.0, 0.0, 0.0;
        sourceNormal1 = source1; // + Eigen::Vector3d::Random();
        sourceNormal1.normalize();
        probe1 << 2.0, 0.0, 0.0;
        probeNormal1 = probe1; // + Eigen::Vector3d::Random();
        probeNormal1.normalize();
        // Reference value
        resultInside = analyticEvaluate(epsInside, sourceNormal1, source1, probeNormal1, probe1);
	    // Evaluation outside the sphere
        source2 << 150.0, 150.0, 150.0;
        sourceNormal2 = source2; // + Eigen::Vector3d::Random();
        sourceNormal2.normalize();
        probe2 << 151.0, 150.0, 150.0;
        probeNormal2 = probe2; // + Eigen::Vector3d::Random();
        probeNormal2.normalize();
        // Reference value
        resultOutside = analyticEvaluate(epsOutside, sourceNormal2, source2, probeNormal2, probe2);
    }
};

BOOST_FIXTURE_TEST_SUITE(TanhSphericalDiffuse2, TanhSphericalDiffuseBogusTest)
/*! \class TanhSphericalDiffuse
 *  \test \b TanhSphericalDiffuseBogusTest_inside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(inside, TanhSphericalDiffuseBogusTest)
{
    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);
    double value = resultInside(0);
    double gf_value = gf.function(source1, probe1);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_CHECK_CLOSE(gf_value, value, 1.0e-07);

    /*
    double derProbe = resultInside(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("derProbe    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_CHECK_CLOSE(derProbe, gf_derProbe, 1.0e-06);

    double derSource = resultInside(2);
    double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
    BOOST_TEST_MESSAGE("derSource    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_CHECK_CLOSE(derSource, gf_derSource, 1.0e-06);
    */
}

/*! \class TanhSphericalDiffuse
 *  \test \b TanhSphericalDiffuseBogusTest_outside tests the evaluation of the TanhSphericalDiffuse Green's function against analytical result for uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(outside, TanhSphericalDiffuseBogusTest)
{
    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius);
    double value = resultOutside(0);
    double gf_value = gf.function(source2, probe2);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("gf_value = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_value);
    BOOST_CHECK_CLOSE(gf_value, value, 1.0e-07);

    /*
    double derProbe = resultOutside(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("derProbe    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derProbe);
    BOOST_TEST_MESSAGE("gf_derProbe = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derProbe);
    BOOST_CHECK_CLOSE(derProbe, gf_derProbe, 1.0e-06);

    double derSource = resultOutside(2);
    double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
    BOOST_TEST_MESSAGE("derSource    = " << std::setprecision(std::numeric_limits<long double>::digits10) << derSource);
    BOOST_TEST_MESSAGE("gf_derSource = " << std::setprecision(std::numeric_limits<long double>::digits10) << gf_derSource);
    BOOST_CHECK_CLOSE(derSource, gf_derSource, 1.0e-06);
    */
}

BOOST_AUTO_TEST_SUITE_END()
