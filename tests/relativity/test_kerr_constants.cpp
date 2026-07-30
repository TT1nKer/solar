#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_constants.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  PASS: " << name << "\n";
        ++passed;
    } else {
        std::cerr << "  FAIL: " << name << "\n";
        ++failed;
    }
}

void check_near(
    const std::string& name,
    double actual,
    double expected,
    double tolerance) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric kerr(1.0, 0.5);
    PhaseSpaceState state{
        0.0,
        Contravariant4{
            Vec4{{0.0, 8.0, 1.0471975511965976, 0.2}}},
        Covariant4{
            Vec4{{-1.0, 0.25, 3.0, 2.0}}},
    };
    check("Kerr constants point is valid", kerr.valid_point(state.x));

    const KerrConstants null_constants =
        evaluate_kerr_constants(
            kerr, state, GeodesicKind::Null);
    check_near(
        "Kerr E from covariant p_t",
        null_constants.E,
        1.0,
        0.0);
    check_near(
        "Kerr Lz from covariant p_phi",
        null_constants.Lz,
        2.0,
        0.0);
    check_near(
        "literal null Carter Q",
        null_constants.Q,
        10.270833333333333,
        2.0e-14);
    check_near(
        "null mass squared",
        null_constants.mass_sq,
        0.0,
        0.0);

    const KerrConstants timelike_constants =
        evaluate_kerr_constants(
            kerr,
            state,
            GeodesicKind::TimelikeUnitMass);
    check_near(
        "literal timelike Carter Q",
        timelike_constants.Q,
        10.333333333333333,
        2.0e-14);
    check_near(
        "timelike mass squared",
        timelike_constants.mass_sq,
        1.0,
        0.0);

    PhaseSpaceState equatorial_zero{};
    equatorial_zero.x = Contravariant4{
        Vec4{{0.0, 8.0, 1.5707963267948966, 0.0}}};
    check_near(
        "equatorial zero-momentum Carter Q",
        evaluate_kerr_constants(
            kerr,
            equatorial_zero,
            GeodesicKind::Null).Q,
        0.0,
        0.0);

    try {
        static_cast<void>(evaluate_kerr_constants(
            kerr,
            state,
            static_cast<GeodesicKind>(91)));
        check("unknown Kerr geodesic kind rejected", false);
    } catch (const std::invalid_argument&) {
        check("unknown Kerr geodesic kind rejected", true);
    }

    PhaseSpaceState non_finite = state;
    non_finite.p.v[2] =
        std::numeric_limits<double>::quiet_NaN();
    try {
        static_cast<void>(evaluate_kerr_constants(
            kerr, non_finite, GeodesicKind::Null));
        check("non-finite Kerr state rejected", false);
    } catch (const std::domain_error&) {
        check("non-finite Kerr state rejected", true);
    }

    PhaseSpaceState invalid_point = state;
    invalid_point.x.v[1] = kerr.outer_horizon_radius();
    try {
        static_cast<void>(evaluate_kerr_constants(
            kerr, invalid_point, GeodesicKind::Null));
        check("invalid Kerr point rejected", false);
    } catch (const std::domain_error&) {
        check("invalid Kerr point rejected", true);
    }

    PhaseSpaceState axis = state;
    axis.x.v[2] = 0.0;
    try {
        static_cast<void>(evaluate_kerr_constants(
            kerr, axis, GeodesicKind::Null));
        check("Kerr polar axis rejected", false);
    } catch (const std::domain_error&) {
        check("Kerr polar axis rejected", true);
    }

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
