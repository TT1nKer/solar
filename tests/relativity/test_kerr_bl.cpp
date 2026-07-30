#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace solar::relativity;

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

void check_relative(const char* name, double actual, double expected,
                    double tolerance) {
    const double scale = std::max(std::fabs(expected), 1.0);
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) / scale <= tolerance);
}

double inverse_identity_error(const Mat4& covariant,
                              const Mat4& contravariant) {
    const Mat4 product = multiply(covariant, contravariant);
    double maximum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            maximum = std::max(
                maximum, std::fabs(product[row][column] - expected));
        }
    }
    return maximum;
}

double matrix_difference(const Mat4& left, const Mat4& right) {
    double maximum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            maximum = std::max(
                maximum, std::fabs(left[row][column] - right[row][column]));
        }
    }
    return maximum;
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric kerr(1.0, 0.9);
    const Contravariant4 equatorial{
        Vec4{{0.0, 10.0, 1.5707963267948966, 0.0}}};

    check("Kerr chart", kerr.chart() == Chart::BoyerLindquist);
    check("Kerr name", kerr.name() == "kerr-bl");
    check_near("Kerr mass", kerr.mass(), 1.0, 0.0);
    check_near("Kerr dimensionless spin", kerr.spin_chi(), 0.9, 0.0);
    check_near("Kerr spin length", kerr.spin_length(), 0.9, 0.0);
    check("ordinary Kerr point accepted", kerr.valid_point(equatorial));

    const Mat4 covariant = kerr.covariant(equatorial);
    const Mat4 contravariant = kerr.contravariant(equatorial);
    check_near("Kerr g_tt", covariant[0][0],
               -0.80000000000000004, 2.0e-15);
    check_near("Kerr g_tphi", covariant[0][3],
               -0.17999999999999999, 2.0e-15);
    check_near("Kerr g_phit", covariant[3][0],
               -0.17999999999999999, 2.0e-15);
    check_near("Kerr g_rr", covariant[1][1],
               1.2374706100730106, 2.0e-15);
    check_near("Kerr g_thetatheta", covariant[2][2],
               100.0, 2.0e-14);
    check_near("Kerr g_phiphi", covariant[3][3],
               100.97200000000001, 2.0e-13);

    check_near("Kerr inverse g_tt", contravariant[0][0],
               -1.2494988244029206, 2.0e-15);
    check_near("Kerr inverse g_tphi", contravariant[0][3],
               -0.0022274470981314192, 2.0e-17);
    check_near("Kerr inverse g_phit", contravariant[3][0],
               -0.0022274470981314192, 2.0e-17);
    check_near("Kerr inverse g_rr", contravariant[1][1],
               0.80810000000000004, 2.0e-15);
    check_near("Kerr inverse g_thetatheta", contravariant[2][2],
               0.01, 2.0e-17);
    check_near("Kerr inverse g_phiphi", contravariant[3][3],
               0.0098997648805840867, 2.0e-17);

    check("Kerr inverse identity below threshold",
          inverse_identity_error(covariant, contravariant) < 5.0e-13);
    check("Kerr covariant symmetry",
          covariant[0][3] == covariant[3][0]);
    check("Kerr inverse symmetry",
          contravariant[0][3] == contravariant[3][0]);
    check("Kerr sampled signature is Lorentzian",
          covariant[1][1] > 0.0 && covariant[2][2] > 0.0 &&
          covariant[0][0] * covariant[3][3] -
              covariant[0][3] * covariant[3][0] < 0.0);

    struct KerrInverseCase {
        double mass;
        double spin_chi;
        Contravariant4 point;
    };
    const KerrInverseCase inverse_cases[] = {
        {1.0, 0.99, Contravariant4{Vec4{{0.0, 2.0, 0.4, 0.0}}}},
        {2.3, -0.7, Contravariant4{Vec4{{1.0, 15.0, 2.0, -0.2}}}},
        {0.4, 0.3, Contravariant4{Vec4{{-2.0, 5.0, 1.3, 0.7}}}},
    };
    double sampled_inverse_error = 0.0;
    bool sampled_signatures_are_lorentzian = true;
    for (const KerrInverseCase& test_case : inverse_cases) {
        const KerrBoyerLindquistMetric sampled_metric(
            test_case.mass, test_case.spin_chi);
        const Mat4 sampled_covariant =
            sampled_metric.covariant(test_case.point);
        const Mat4 sampled_contravariant =
            sampled_metric.contravariant(test_case.point);
        sampled_inverse_error = std::max(
            sampled_inverse_error,
            inverse_identity_error(
                sampled_covariant, sampled_contravariant));
        sampled_signatures_are_lorentzian &=
            sampled_covariant[1][1] > 0.0 &&
            sampled_covariant[2][2] > 0.0 &&
            sampled_covariant[0][0] * sampled_covariant[3][3] -
                sampled_covariant[0][3] * sampled_covariant[3][0] < 0.0;
    }
    check("Kerr inverse identity across sampled exterior points",
          sampled_inverse_error < 5.0e-13);
    check("Kerr signature across sampled exterior points",
          sampled_signatures_are_lorentzian);
    std::cout << "    sampled max inverse error: "
              << sampled_inverse_error << '\n';

    check_near("Kerr outer horizon", kerr.outer_horizon_radius(),
               1.4358898943540672, 1.0e-15);
    check_near("Kerr inner horizon", kerr.inner_horizon_radius(),
               0.5641101056459328, 1.0e-15);
    check_near("Kerr equatorial stationary limit",
               kerr.outer_stationary_limit_radius(
                   1.5707963267948966),
               2.0, 1.0e-15);
    check_near("Kerr polar stationary limit",
               kerr.outer_stationary_limit_radius(0.0),
               kerr.outer_horizon_radius(), 1.0e-15);
    try {
        (void)kerr.outer_stationary_limit_radius(
            std::numeric_limits<double>::quiet_NaN());
        check("non-finite stationary-limit angle rejected", false);
    } catch (const std::invalid_argument&) {
        check("non-finite stationary-limit angle rejected", true);
    }

    const double large_mass = 1.0e200;
    const KerrBoyerLindquistMetric large_kerr(large_mass, 0.5);
    check_relative(
        "large finite mass has stable outer horizon",
        large_kerr.outer_horizon_radius(),
        large_mass * (1.0 + std::sqrt(0.75)),
        2.0e-15);
    check_relative(
        "large finite mass has stable inner horizon",
        large_kerr.inner_horizon_radius(),
        large_mass * (1.0 - std::sqrt(0.75)),
        2.0e-15);
    const double stationary_theta = 0.7;
    check_relative(
        "large finite mass has stable stationary limit",
        large_kerr.outer_stationary_limit_radius(stationary_theta),
        large_mass * (
            1.0 + std::sqrt(
                1.0 -
                0.25 * std::cos(stationary_theta) *
                    std::cos(stationary_theta))),
        2.0e-15);

    const KerrBoyerLindquistMetric zero_spin(1.7, 0.0);
    const SchwarzschildBoyerLindquistMetric schwarzschild(1.7);
    for (const Contravariant4 point : {
             Contravariant4{Vec4{{0.0, 6.0, 0.7, -1.0}}},
             Contravariant4{Vec4{{4.0, 20.0, 1.2, 2.0}}},
         }) {
        check("Kerr a=0 covariant equals Schwarzschild",
              matrix_difference(
                  zero_spin.covariant(point),
                  schwarzschild.covariant(point)) < 2.0e-13);
        check("Kerr a=0 inverse equals Schwarzschild",
              matrix_difference(
                  zero_spin.contravariant(point),
                  schwarzschild.contravariant(point)) < 2.0e-15);
    }

    Contravariant4 horizon = equatorial;
    horizon.v[1] = kerr.outer_horizon_radius();
    check("Kerr horizon rejected", !kerr.valid_point(horizon));
    Contravariant4 axis = equatorial;
    axis.v[2] = 0.0;
    check("Kerr axis rejected", !kerr.valid_point(axis));
    Contravariant4 non_finite = equatorial;
    non_finite.v[1] = std::numeric_limits<double>::infinity();
    check("non-finite Kerr point rejected", !kerr.valid_point(non_finite));
    Contravariant4 overflowing_A = equatorial;
    overflowing_A.v[1] = 1.0e100;
    check("Kerr radius with overflowing A rejected",
          !kerr.valid_point(overflowing_A));

    try {
        (void)KerrBoyerLindquistMetric(0.0, 0.5);
        check("non-positive Kerr mass rejected", false);
    } catch (const std::invalid_argument&) {
        check("non-positive Kerr mass rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(
            std::numeric_limits<double>::quiet_NaN(), 0.5);
        check("non-finite Kerr mass rejected", false);
    } catch (const std::invalid_argument&) {
        check("non-finite Kerr mass rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(
            1.0, std::numeric_limits<double>::quiet_NaN());
        check("non-finite Kerr spin rejected", false);
    } catch (const std::invalid_argument&) {
        check("non-finite Kerr spin rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(1.0, 0.5, 0.0);
        check("non-positive Kerr margin rejected", false);
    } catch (const std::invalid_argument&) {
        check("non-positive Kerr margin rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(
            std::numeric_limits<double>::denorm_min(), 0.5);
        check("underflowing Kerr margin rejected", false);
    } catch (const std::invalid_argument&) {
        check("underflowing Kerr margin rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(1.0, 1.0);
        check("extremal Kerr spin rejected", false);
    } catch (const std::invalid_argument&) {
        check("extremal Kerr spin rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(1.0e308, 0.0);
        check("unrepresentable Kerr horizon rejected", false);
    } catch (const std::invalid_argument&) {
        check("unrepresentable Kerr horizon rejected", true);
    }
    try {
        (void)KerrBoyerLindquistMetric(1.0e200, 0.5, 1.0e200);
        check("unrepresentable Kerr margin rejected", false);
    } catch (const std::invalid_argument&) {
        check("unrepresentable Kerr margin rejected", true);
    }
    try {
        (void)kerr.covariant(horizon);
        check("invalid Kerr evaluation throws", false);
    } catch (const std::domain_error&) {
        check("invalid Kerr evaluation throws", true);
    }

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
