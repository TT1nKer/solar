#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/schwarzschild_metric.h"

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

double derivative_max_abs(const std::array<Mat4, 4>& derivatives) {
    double maximum = 0.0;
    for (const auto& derivative : derivatives) {
        for (const auto& row : derivative) {
            for (const double component : row) {
                maximum = std::max(maximum, std::fabs(component));
            }
        }
    }
    return maximum;
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

} // namespace

int main() {
    const MinkowskiMetric metric;
    const Contravariant4 x{Vec4{{2.0, -3.0, 4.0, 5.0}}};

    check("Minkowski chart", metric.chart() == Chart::MinkowskiCartesian);
    check("Minkowski name", metric.name() == "minkowski");
    check("finite point accepted", metric.valid_point(x));

    const Mat4 covariant = metric.covariant(x);
    const Mat4 contravariant = metric.contravariant(x);
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected =
                row == column ? (row == 0 ? -1.0 : 1.0) : 0.0;
            check_near("Minkowski covariant component",
                       covariant[row][column], expected, 0.0);
            check_near("Minkowski inverse component",
                       contravariant[row][column], expected, 0.0);
        }
    }

    check_near("Minkowski derivatives are exactly zero",
               derivative_max_abs(metric.contravariant_derivatives(x)),
               0.0, 0.0);

    Contravariant4 invalid = x;
    invalid.v[2] = std::numeric_limits<double>::quiet_NaN();
    check("non-finite Minkowski point rejected", !metric.valid_point(invalid));
    try {
        (void)metric.covariant(invalid);
        check("invalid Minkowski evaluation throws", false);
    } catch (const std::domain_error&) {
        check("invalid Minkowski evaluation throws", true);
    }

    const SchwarzschildBoyerLindquistMetric schwarzschild(1.0);
    const Contravariant4 exterior{
        Vec4{{0.0, 10.0, 1.5707963267948966, 0.0}}};
    check("Schwarzschild chart",
          schwarzschild.chart() == Chart::BoyerLindquist);
    check("Schwarzschild name",
          schwarzschild.name() == "schwarzschild");
    check_near("Schwarzschild horizon", schwarzschild.outer_horizon_radius(),
               2.0, 0.0);
    check("ordinary Schwarzschild point accepted",
          schwarzschild.valid_point(exterior));

    const Mat4 schwarzschild_covariant = schwarzschild.covariant(exterior);
    const Mat4 schwarzschild_inverse =
        schwarzschild.contravariant(exterior);
    check_near("Schwarzschild g_tt",
               schwarzschild_covariant[0][0], -0.8, 1.0e-15);
    check_near("Schwarzschild g_rr",
               schwarzschild_covariant[1][1], 1.25, 1.0e-15);
    check_near("Schwarzschild g_thetatheta",
               schwarzschild_covariant[2][2], 100.0, 1.0e-14);
    check_near("Schwarzschild g_phiphi",
               schwarzschild_covariant[3][3], 100.0, 1.0e-14);
    check_near("Schwarzschild inverse g_tt",
               schwarzschild_inverse[0][0], -1.25, 1.0e-15);
    check_near("Schwarzschild inverse g_rr",
               schwarzschild_inverse[1][1], 0.8, 1.0e-15);
    check_near("Schwarzschild inverse g_thetatheta",
               schwarzschild_inverse[2][2], 0.01, 1.0e-17);
    check_near("Schwarzschild inverse g_phiphi",
               schwarzschild_inverse[3][3], 0.01, 1.0e-17);
    check("Schwarzschild inverse identity",
          inverse_identity_error(
              schwarzschild_covariant, schwarzschild_inverse) < 5.0e-13);

    const auto schwarzschild_covariant_derivatives =
        schwarzschild.covariant_derivatives(exterior);
    check_near("Schwarzschild dr g_tt",
               schwarzschild_covariant_derivatives[1][0][0],
               -0.02, 1.0e-16);
    check_near("Schwarzschild dr g_rr",
               schwarzschild_covariant_derivatives[1][1][1],
               -0.03125, 1.0e-16);
    check_near("Schwarzschild dr g_thetatheta",
               schwarzschild_covariant_derivatives[1][2][2],
               20.0, 1.0e-14);
    check_near("Schwarzschild dr g_phiphi",
               schwarzschild_covariant_derivatives[1][3][3],
               20.0, 1.0e-14);

    const auto schwarzschild_inverse_derivatives =
        schwarzschild.contravariant_derivatives(exterior);
    check_near("Schwarzschild dr inverse g_tt",
               schwarzschild_inverse_derivatives[1][0][0],
               0.03125, 1.0e-15);
    check_near("Schwarzschild dr inverse g_rr",
               schwarzschild_inverse_derivatives[1][1][1],
               0.02, 1.0e-16);
    check_near("Schwarzschild dr inverse g_thetatheta",
               schwarzschild_inverse_derivatives[1][2][2],
               -0.002, 1.0e-17);
    check_near("Schwarzschild dr inverse g_phiphi",
               schwarzschild_inverse_derivatives[1][3][3],
               -0.002, 1.0e-17);

    Contravariant4 horizon = exterior;
    horizon.v[1] = 2.0;
    check("Schwarzschild horizon rejected",
          !schwarzschild.valid_point(horizon));
    Contravariant4 axis = exterior;
    axis.v[2] = 0.0;
    check("Schwarzschild polar axis rejected",
          !schwarzschild.valid_point(axis));
    Contravariant4 outside_theta = exterior;
    outside_theta.v[2] = 4.0;
    check("Schwarzschild theta outside chart rejected",
          !schwarzschild.valid_point(outside_theta));
    Contravariant4 unrepresentable_radius = exterior;
    unrepresentable_radius.v[1] = 1.0e200;
    check("Schwarzschild radius with overflowing r squared rejected",
          !schwarzschild.valid_point(unrepresentable_radius));
    try {
        (void)schwarzschild.covariant(unrepresentable_radius);
        check("unrepresentable Schwarzschild evaluation throws", false);
    } catch (const std::domain_error&) {
        check("unrepresentable Schwarzschild evaluation throws", true);
    }

    try {
        (void)SchwarzschildBoyerLindquistMetric(0.0);
        check("non-positive Schwarzschild mass rejected", false);
    } catch (const std::invalid_argument&) {
        check("non-positive Schwarzschild mass rejected", true);
    }
    try {
        (void)SchwarzschildBoyerLindquistMetric(1.0e308);
        check("unrepresentable Schwarzschild horizon rejected", false);
    } catch (const std::invalid_argument&) {
        check("unrepresentable Schwarzschild horizon rejected", true);
    }
    try {
        (void)SchwarzschildBoyerLindquistMetric(
            1.0e200, 1.0e200);
        check("unrepresentable Schwarzschild margin rejected", false);
    } catch (const std::invalid_argument&) {
        check("unrepresentable Schwarzschild margin rejected", true);
    }

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
