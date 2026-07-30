#include "solar/relativity/kerr_constants.h"

#include "../../src/relativity/kerr_separated_potentials.h"

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

template <typename Function>
void check_invalid_argument(
    const std::string& name,
    Function&& function) {
    try {
        function();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    }
}

template <typename Function>
void check_domain_error(
    const std::string& name,
    Function&& function) {
    try {
        function();
        check(name, false);
    } catch (const std::domain_error&) {
        check(name, true);
    }
}

} // namespace

int main() {
    const KerrConstants null_constants{1.0, 2.0, 3.0, 0.0};
    const detail::KerrSeparatedPotentials null_potential(
        1.0, 0.5, null_constants);
    const auto values = null_potential.evaluate(4.0, 0.25);

    check_near("literal Delta", values.delta, 8.25, 1.0e-15);
    check_near(
        "literal Sigma", values.sigma, 16.015625, 1.0e-15);
    check_near(
        "literal radial R", values.radial, 189.25, 1.0e-13);
    check_near(
        "literal polar U",
        values.polar,
        2.5771484375,
        1.0e-14);
    check_near(
        "literal radial derivative",
        values.radial_derivative,
        212.5,
        1.0e-13);
    check_near(
        "literal polar derivative",
        values.polar_derivative,
        -3.390625,
        1.0e-14);
    check_near(
        "radial scale includes potential",
        values.radial_scale,
        189.25,
        1.0e-13);
    check_near(
        "polar scale includes potential",
        values.polar_scale,
        2.5771484375,
        1.0e-14);

    const detail::KerrSeparatedPotentials timelike_potential(
        1.0,
        0.5,
        KerrConstants{1.0, 2.0, 3.0, 1.0});
    const auto timelike = timelike_potential.evaluate(4.0, 0.25);
    check_near(
        "literal timelike radial R",
        timelike.radial,
        57.25,
        1.0e-13);
    check_near(
        "literal timelike radial derivative",
        timelike.radial_derivative,
        50.5,
        1.0e-13);
    check_near(
        "literal timelike polar U",
        timelike.polar,
        2.5625,
        1.0e-14);
    check_near(
        "literal timelike polar derivative",
        timelike.polar_derivative,
        -3.5,
        1.0e-14);

    const detail::KerrSeparatedPotentials negative_spin(
        1.0, -0.5, null_constants);
    const auto negative_spin_values =
        negative_spin.evaluate(4.0, 0.25);
    check_near(
        "negative spin changes radial potential with fixed Lz",
        negative_spin_values.radial,
        221.25,
        1.0e-13);
    check_near(
        "negative spin radial derivative",
        negative_spin_values.radial_derivative,
        220.5,
        1.0e-13);
    check_near(
        "polar potential depends on spin squared",
        negative_spin_values.polar,
        values.polar,
        0.0);

    const double radial_h = 1.0e-4;
    const double radial_finite_difference =
        (-null_potential.evaluate(4.0 + 2.0 * radial_h, 0.25)
              .radial +
         8.0 *
             null_potential.evaluate(4.0 + radial_h, 0.25)
                 .radial -
         8.0 *
             null_potential.evaluate(4.0 - radial_h, 0.25)
                 .radial +
         null_potential.evaluate(4.0 - 2.0 * radial_h, 0.25)
             .radial) /
        (12.0 * radial_h);
    check(
        "radial derivative matches five-point difference",
        std::fabs(
            radial_finite_difference -
            values.radial_derivative) /
                std::fabs(values.radial_derivative) <
            1.0e-10);

    const double polar_h = 1.0e-4;
    const double polar_finite_difference =
        (-null_potential.evaluate(4.0, 0.25 + 2.0 * polar_h)
              .polar +
         8.0 *
             null_potential.evaluate(4.0, 0.25 + polar_h)
                 .polar -
         8.0 *
             null_potential.evaluate(4.0, 0.25 - polar_h)
                 .polar +
         null_potential.evaluate(4.0, 0.25 - 2.0 * polar_h)
             .polar) /
        (12.0 * polar_h);
    check(
        "polar derivative matches five-point difference",
        std::fabs(
            polar_finite_difference -
            values.polar_derivative) /
                std::fabs(values.polar_derivative) <
            1.0e-10);

    check_invalid_argument(
        "non-positive mass rejected",
        [&] {
            static_cast<void>(
                detail::KerrSeparatedPotentials(
                    0.0, 0.5, null_constants));
        });
    check_invalid_argument(
        "non-finite spin rejected",
        [&] {
            static_cast<void>(
                detail::KerrSeparatedPotentials(
                    1.0,
                    std::numeric_limits<double>::infinity(),
                    null_constants));
        });
    check_invalid_argument(
        "non-finite constants rejected",
        [&] {
            KerrConstants invalid = null_constants;
            invalid.Q =
                std::numeric_limits<double>::quiet_NaN();
            static_cast<void>(
                detail::KerrSeparatedPotentials(
                    1.0, 0.5, invalid));
        });
    check_domain_error(
        "mu outside cosine domain rejected",
        [&] {
            static_cast<void>(
                null_potential.evaluate(4.0, 1.01));
        });
    check_domain_error(
        "non-finite radius rejected",
        [&] {
            static_cast<void>(null_potential.evaluate(
                std::numeric_limits<double>::infinity(),
                0.25));
        });
    check_domain_error(
        "overflowing potential rejected",
        [&] {
            static_cast<void>(
                null_potential.evaluate(1.0e200, 0.25));
        });

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
