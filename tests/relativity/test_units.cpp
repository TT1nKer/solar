#include "solar/relativity/units.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using solar::relativity::GeometricUnits;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(const char* name, double actual, double expected,
                double tolerance) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

template <typename Operation>
void check_invalid_argument(const char* name, Operation operation) {
    try {
        operation();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    } catch (...) {
        check(name, false);
    }
}

} // namespace

int main() {
    const auto units = GeometricUnits::from_solar_masses(1.0);

    check_near("solar mass in kilograms", units.mass_kg, 1.98847e30, 1.0e20);
    check_near("one solar mass length scale", units.M_length_m,
               1476.6696910334392, 1.0e-9);
    check_near("one solar mass time scale", units.M_time_s,
               4.925639893961039e-6, 1.0e-18);

    check_near("length SI to M", units.length_si_to_M(units.M_length_m * 3.5),
               3.5, 1.0e-14);
    check_near("length round trip",
               units.length_M_to_si(units.length_si_to_M(12345.0)),
               12345.0, 1.0e-10);
    check_near("time SI to M", units.time_si_to_M(units.M_time_s * 7.25),
               7.25, 1.0e-14);
    check_near("time round trip",
               units.time_M_to_si(units.time_si_to_M(0.125)),
               0.125, 1.0e-15);
    check_near("speed of light is one", units.velocity_si_to_c(299792458.0),
               1.0, 1.0e-15);

    check_invalid_argument("zero mass rejected", [] {
        (void)GeometricUnits::from_mass_kg(0.0);
    });
    check_invalid_argument("negative solar masses rejected", [] {
        (void)GeometricUnits::from_solar_masses(-1.0);
    });
    check_invalid_argument("non-finite mass rejected", [] {
        (void)GeometricUnits::from_mass_kg(
            std::numeric_limits<double>::infinity());
    });
    check_invalid_argument("non-finite conversion rejected", [&units] {
        (void)units.length_si_to_M(
            std::numeric_limits<double>::quiet_NaN());
    });
    check_invalid_argument("uninitialized length scale rejected", [] {
        (void)GeometricUnits{}.length_si_to_M(1.0);
    });
    check_invalid_argument("uninitialized time scale rejected", [] {
        (void)GeometricUnits{}.time_M_to_si(1.0);
    });

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
