#include "solar/relativity/kerr_orbits.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <iomanip>
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
    check(
        name,
        std::isfinite(actual) &&
            std::fabs(actual - expected) <= tolerance);
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric schwarzschild(2.0, 0.0);
    for (const OrbitSense sense :
         {OrbitSense::Prograde, OrbitSense::Retrograde}) {
        check_near(
            "Schwarzschild ISCO radius",
            kerr_isco_radius(schwarzschild, sense),
            12.0,
            0.0);
        check_near(
            "Schwarzschild photon radius",
            kerr_equatorial_photon_radius(
                schwarzschild, sense),
            6.0,
            0.0);
        check_near(
            "Schwarzschild marginally bound radius",
            kerr_marginally_bound_radius(
                schwarzschild, sense),
            8.0,
            0.0);
    }

    const KerrBoyerLindquistMetric positive_spin(1.0, 0.5);
    check_near(
        "positive-spin prograde ISCO",
        kerr_isco_radius(
            positive_spin, OrbitSense::Prograde),
        4.233002529530826,
        2.0e-15);
    check_near(
        "positive-spin retrograde ISCO",
        kerr_isco_radius(
            positive_spin, OrbitSense::Retrograde),
        7.554584714512358,
        3.0e-15);
    check_near(
        "positive-spin prograde photon radius",
        kerr_equatorial_photon_radius(
            positive_spin, OrbitSense::Prograde),
        2.347296355333861,
        2.0e-15);
    check_near(
        "positive-spin retrograde photon radius",
        kerr_equatorial_photon_radius(
            positive_spin, OrbitSense::Retrograde),
        3.532088886237956,
        2.0e-15);
    check_near(
        "positive-spin prograde marginally bound",
        kerr_marginally_bound_radius(
            positive_spin, OrbitSense::Prograde),
        2.914213562373095,
        2.0e-15);
    check_near(
        "positive-spin retrograde marginally bound",
        kerr_marginally_bound_radius(
            positive_spin, OrbitSense::Retrograde),
        4.949489742783178,
        2.0e-15);

    const KerrBoyerLindquistMetric negative_spin(1.0, -0.5);
    check_near(
        "negative-spin prograde radius remains relative to spin",
        kerr_isco_radius(
            negative_spin, OrbitSense::Prograde),
        4.233002529530826,
        2.0e-15);
    check_near(
        "negative-spin retrograde radius remains relative to spin",
        kerr_isco_radius(
            negative_spin, OrbitSense::Retrograde),
        7.554584714512358,
        3.0e-15);

    try {
        static_cast<void>(kerr_isco_radius(
            positive_spin, static_cast<OrbitSense>(27)));
        check("unknown orbit sense rejected", false);
    } catch (const std::invalid_argument&) {
        check("unknown orbit sense rejected", true);
    }

    const KerrBoyerLindquistMetric huge_kerr(
        std::numeric_limits<double>::max() / 2.0,
        0.5);
    try {
        static_cast<void>(kerr_isco_radius(
            huge_kerr, OrbitSense::Retrograde));
        check("overflowing Kerr ISCO rejected", false);
    } catch (const std::overflow_error&) {
        check("overflowing Kerr ISCO rejected", true);
    }
    try {
        static_cast<void>(kerr_equatorial_photon_radius(
            huge_kerr, OrbitSense::Retrograde));
        check("overflowing Kerr photon radius rejected", false);
    } catch (const std::overflow_error&) {
        check("overflowing Kerr photon radius rejected", true);
    }
    try {
        static_cast<void>(kerr_marginally_bound_radius(
            huge_kerr, OrbitSense::Retrograde));
        check("overflowing Kerr marginal radius rejected", false);
    } catch (const std::overflow_error&) {
        check("overflowing Kerr marginal radius rejected", true);
    }

    const KerrBoyerLindquistMetric unit_schwarzschild(
        1.0, 0.0);
    const CircularOrbitResult prograde_six =
        evaluate_equatorial_circular_timelike_orbit(
            unit_schwarzschild, 6.0, OrbitSense::Prograde);
    check("Schwarzschild r=6 circular orbit exists",
          bool(prograde_six));
    check_near(
        "Schwarzschild circular specific energy",
        prograde_six.orbit->specific_energy,
        0.9428090415820634,
        2.0e-15);
    check_near(
        "Schwarzschild prograde specific Lz",
        prograde_six.orbit->specific_lz,
        3.4641016151377544,
        3.0e-15);
    check_near(
        "Schwarzschild circular angular velocity",
        prograde_six.orbit->angular_velocity,
        0.06804138174397717,
        2.0e-17);
    check(
        "Schwarzschild ISCO is classified stable",
        prograde_six.orbit->stability ==
            CircularOrbitStability::Stable);

    const CircularOrbitResult retrograde_six =
        evaluate_equatorial_circular_timelike_orbit(
            unit_schwarzschild, 6.0, OrbitSense::Retrograde);
    check_near(
        "Schwarzschild retrograde specific Lz",
        retrograde_six.orbit->specific_lz,
        -3.4641016151377544,
        3.0e-15);
    check_near(
        "Schwarzschild retrograde angular velocity",
        retrograde_six.orbit->angular_velocity,
        -0.06804138174397717,
        2.0e-17);

    const CircularOrbitResult unstable =
        evaluate_equatorial_circular_timelike_orbit(
            unit_schwarzschild, 5.0, OrbitSense::Prograde);
    check("Schwarzschild r=5 timelike circular orbit exists",
          bool(unstable));
    check(
        "Schwarzschild r=5 is unstable",
        unstable.orbit->stability ==
            CircularOrbitStability::Unstable);
    check(
        "Schwarzschild photon radius is not timelike",
        !evaluate_equatorial_circular_timelike_orbit(
            unit_schwarzschild,
            3.0,
            OrbitSense::Prograde));

    const CircularOrbitResult negative_prograde =
        evaluate_equatorial_circular_timelike_orbit(
            negative_spin, 8.0, OrbitSense::Prograde);
    const CircularOrbitResult negative_retrograde =
        evaluate_equatorial_circular_timelike_orbit(
            negative_spin, 8.0, OrbitSense::Retrograde);
    check("negative-spin prograde orbit exists",
          bool(negative_prograde));
    check("negative-spin retrograde orbit exists",
          bool(negative_retrograde));
    check("negative-spin prograde Omega is negative",
          negative_prograde.orbit->angular_velocity < 0.0);
    check("negative-spin prograde Lz is negative",
          negative_prograde.orbit->specific_lz < 0.0);
    check("negative-spin retrograde Omega is positive",
          negative_retrograde.orbit->angular_velocity > 0.0);
    check("negative-spin retrograde Lz is positive",
          negative_retrograde.orbit->specific_lz > 0.0);

    const Contravariant4 circular_point{
        Vec4{{0.0, 8.0, 1.5707963267948966, 0.0}}};
    const ObserverResult circular_observer =
        make_equatorial_circular_observer(
            positive_spin,
            circular_point,
            OrbitSense::Prograde);
    check("Kerr circular observer exists", bool(circular_observer));
    const CircularOrbitResult circular_properties =
        evaluate_equatorial_circular_timelike_orbit(
            positive_spin, 8.0, OrbitSense::Prograde);
    const Covariant4 lowered_velocity = lower_index(
        positive_spin.covariant(circular_point),
        circular_observer.frame->tetrad.basis[0]);
    check_near(
        "circular energy matches lowered observer",
        -lowered_velocity.v[0],
        circular_properties.orbit->specific_energy,
        2.0e-13);
    check_near(
        "circular Lz matches lowered observer",
        lowered_velocity.v[3],
        circular_properties.orbit->specific_lz,
        2.0e-13);
    check(
        "circular observer tetrad gate",
        tetrad_orthonormality_error(
            positive_spin, *circular_observer.frame) < 1.0e-10);

    Contravariant4 non_equatorial = circular_point;
    non_equatorial.v[2] = 1.2;
    check(
        "non-equatorial circular observer rejected",
        make_equatorial_circular_observer(
            positive_spin,
            non_equatorial,
            OrbitSense::Prograde).error ==
            ObserverError::CircularWorldlineNotTimelike);

    Contravariant4 photon_radius_point = circular_point;
    photon_radius_point.v[1] =
        kerr_equatorial_photon_radius(
            positive_spin, OrbitSense::Prograde);
    check(
        "photon-radius circular observer rejected",
        make_equatorial_circular_observer(
            positive_spin,
            photon_radius_point,
            OrbitSense::Prograde).error ==
            ObserverError::CircularWorldlineNotTimelike);

    std::cout << std::setprecision(17)
              << "  prograde_isco="
              << kerr_isco_radius(
                     positive_spin, OrbitSense::Prograde)
              << " retrograde_isco="
              << kerr_isco_radius(
                     positive_spin, OrbitSense::Retrograde)
              << " prograde_photon="
              << kerr_equatorial_photon_radius(
                     positive_spin, OrbitSense::Prograde)
              << " retrograde_photon="
              << kerr_equatorial_photon_radius(
                     positive_spin, OrbitSense::Retrograde)
              << " circular_energy_error="
              << std::fabs(
                     -lowered_velocity.v[0] -
                     circular_properties.orbit->specific_energy)
              << " circular_lz_error="
              << std::fabs(
                     lowered_velocity.v[3] -
                     circular_properties.orbit->specific_lz)
              << " circular_tetrad_error="
              << tetrad_orthonormality_error(
                     positive_spin, *circular_observer.frame)
              << "\n";
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
