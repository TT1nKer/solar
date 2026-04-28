#include "solar/kepler.h"
#include "solar/ephemeris.h"
#include "solar/transfer.h"
#include "solar/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace solar;
using namespace solar::constants;

static int tests_passed = 0;
static int tests_failed = 0;

static void check(const char* name, double got, double expected, double tol) {
    double err = std::fabs(got - expected);
    bool ok = err < tol;
    if (ok) {
        std::cout << "  PASS: " << name << " = " << got << "\n";
        tests_passed++;
    } else {
        std::cout << "  FAIL: " << name << " = " << got
                  << " (expected " << expected << ", err " << err << ")\n";
        tests_failed++;
    }
}

static void test_kepler_solver() {
    std::cout << "=== Kepler Equation Solver ===\n";

    // E - 0.5*sin(E) = 1.0 => E ~ 1.4987 (verified numerically)
    double E = solve_kepler(1.0, 0.5);
    check("solve_kepler(M=1.0, e=0.5)", E, 1.49870113, 1e-6);

    // Circular orbit: e=0, E should equal M
    E = solve_kepler(2.0, 0.0);
    check("solve_kepler(M=2.0, e=0.0)", E, 2.0, 1e-10);

    // High eccentricity
    E = solve_kepler(1.0, 0.9);
    check("solve_kepler(M=1.0, e=0.9)", E, 1.862087, 1e-4);
}

static void test_orbital_period() {
    std::cout << "\n=== Orbital Period ===\n";

    // Earth: a = 1 AU, period should be ~365.25 days
    double a_earth = 1.00000261 * AU;
    double T = orbital_period(a_earth, MU_SUN);
    double T_days = T / DAY;
    check("Earth period (days)", T_days, 365.26, 0.1);

    // Mars: period ~687 days
    double a_mars = 1.52371034 * AU;
    T = orbital_period(a_mars, MU_SUN);
    T_days = T / DAY;
    check("Mars period (days)", T_days, 687.0, 1.0);
}

static void test_julian_date() {
    std::cout << "\n=== Julian Date Conversion ===\n";

    // J2000.0 = 2000-01-01 12:00 = JD 2451545.0
    double jd = calendar_to_jd(2000, 1, 1, 12.0);
    check("J2000.0 JD", jd, 2451545.0, 0.001);

    // Round-trip test
    int y, m, d;
    double h;
    jd_to_calendar(2451545.0, y, m, d, h);
    check("J2000 year", y, 2000, 0.5);
    check("J2000 month", m, 1, 0.5);
    check("J2000 day", d, 1, 0.5);
}

static void test_elements_roundtrip() {
    std::cout << "\n=== Elements <-> State Roundtrip ===\n";

    // Create Earth-like elements
    OrbitalElements elem;
    elem.a = 1.0 * AU;
    elem.e = 0.0167;
    elem.i = 0.0 * DEG2RAD;
    elem.omega = 102.9 * DEG2RAD;
    elem.Omega = 0.0;
    elem.M0 = 100.0 * DEG2RAD;
    elem.epoch = J2000;

    double nu = mean_to_true_anomaly(elem.M0, elem.e);
    State s = elements_to_state(elem, MU_SUN, nu);

    // Convert back
    OrbitalElements elem2 = state_to_elements(s, MU_SUN, J2000);

    check("roundtrip a", elem2.a, elem.a, elem.a * 1e-8);
    check("roundtrip e", elem2.e, elem.e, 1e-8);
}

static void test_hohmann() {
    std::cout << "\n=== Hohmann Transfer (Earth -> Mars) ===\n";

    double r1 = 1.0 * AU;    // Earth orbit radius
    double r2 = 1.524 * AU;  // Mars orbit radius

    TransferResult result = hohmann_transfer(r1, r2, MU_SUN);

    check("dv_total (km/s)", result.dv_total, 5.59, 0.1);
    check("TOF (days)", result.tof / DAY, 259.0, 5.0);
}

static void test_ephemeris() {
    std::cout << "\n=== Ephemeris ===\n";

    auto bodies = load_solar_system();

    // Check that we loaded all bodies (9 planets+sun + 8 moons = 17)
    check("body count", bodies.size(), 17, 0.5);

    // Earth at J2000 should be ~1 AU from Sun
    const Body* earth = find_body(bodies, "Earth");
    if (earth) {
        State s = compute_position(*earth, J2000);
        double dist_au = s.pos.norm() / AU;
        check("Earth distance at J2000 (AU)", dist_au, 1.0, 0.02);

        double speed = s.vel.norm();
        check("Earth speed at J2000 (km/s)", speed, 29.78, 1.0);
    } else {
        std::cout << "  FAIL: Earth not found\n";
        tests_failed++;
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    test_kepler_solver();
    test_orbital_period();
    test_julian_date();
    test_elements_roundtrip();
    test_hohmann();
    test_ephemeris();

    std::cout << "\n=== Results: " << tests_passed << " passed, "
              << tests_failed << " failed ===\n";

    return tests_failed > 0 ? 1 : 0;
}
