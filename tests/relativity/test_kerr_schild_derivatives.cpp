#include "solar/relativity/dual4.h"
#include "solar/relativity/kerr_schild_metric.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(
    const char* name,
    bool condition,
    double diagnostic = 0.0) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name
              << " (max normalized error=" << diagnostic << ")\n";
    condition ? ++passed : ++failed;
}

double normalized_difference(
    double actual,
    double expected) {
    return std::fabs(actual - expected) /
           std::max(
               {1.0, std::fabs(actual), std::fabs(expected)});
}

Contravariant4 from_oblate(
    double radius,
    double theta,
    double azimuth,
    double spin_a) {
    const double sine = std::sin(theta);
    return Contravariant4{Vec4{{
        0.7,
        (radius * std::cos(azimuth) -
         spin_a * std::sin(azimuth)) * sine,
        (radius * std::sin(azimuth) +
         spin_a * std::cos(azimuth)) * sine,
        radius * std::cos(theta),
    }}};
}

Vec3 independent_dual_radius_gradient(
    const Contravariant4& point,
    double spin_a) {
    const Dual4 x = Dual4::variable(point.v[1], 1);
    const Dual4 y = Dual4::variable(point.v[2], 2);
    const Dual4 z = Dual4::variable(point.v[3], 3);
    const double spin_squared = spin_a * spin_a;
    const Dual4 rho_squared = x * x + y * y + z * z;
    const Dual4 q = rho_squared - spin_squared;
    const Dual4 radius_squared =
        0.5 * (
            q + sqrt(
                q * q +
                4.0 * spin_squared * z * z));
    const Dual4 radius = sqrt(radius_squared);
    return Vec3{{
        radius.derivative[1],
        radius.derivative[2],
        radius.derivative[3],
    }};
}

Mat4 five_point_inverse_derivative(
    const KerrSchildCartesianMetric& metric,
    const Contravariant4& point,
    std::size_t coordinate,
    double step) {
    Contravariant4 plus_two = point;
    Contravariant4 plus_one = point;
    Contravariant4 minus_one = point;
    Contravariant4 minus_two = point;
    plus_two.v[coordinate] += 2.0 * step;
    plus_one.v[coordinate] += step;
    minus_one.v[coordinate] -= step;
    minus_two.v[coordinate] -= 2.0 * step;

    const Mat4 f_plus_two = metric.contravariant(plus_two);
    const Mat4 f_plus_one = metric.contravariant(plus_one);
    const Mat4 f_minus_one = metric.contravariant(minus_one);
    const Mat4 f_minus_two = metric.contravariant(minus_two);
    Mat4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result[row][column] =
                (-f_plus_two[row][column] +
                 8.0 * f_plus_one[row][column] -
                 8.0 * f_minus_one[row][column] +
                 f_minus_two[row][column]) /
                (12.0 * step);
        }
    }
    return result;
}

struct DerivativeFixture {
    KerrSchildCartesianMetric metric;
    Contravariant4 point;
};

} // namespace

int main() {
    const KerrSchildCartesianMetric positive(1.0, 0.8);
    const KerrSchildCartesianMetric negative(2.0, -0.65);
    const std::vector<DerivativeFixture> fixtures{
        {positive, from_oblate(8.0, 1.1, 0.4, 0.8)},
        {positive,
         Contravariant4{
             Vec4{{-1.0, 1.0e-8, -2.0e-8, 2.0}}}},
        {positive,
         from_oblate(
             positive.outer_horizon_radius() + 1.0e-4,
             1.3,
             -0.7,
             0.8)},
        {positive,
         Contravariant4{
             Vec4{{2.0, 0.2, 0.1, 0.3}}}},
        {negative, from_oblate(9.0, 2.0, 1.2, -1.3)},
        {negative, from_oblate(1.0, 0.8, -0.5, -1.3)},
    };

    double maximum_gradient_error = 0.0;
    double maximum_derivative_error = 0.0;
    std::size_t worst_fixture = 0;
    std::size_t worst_coordinate = 0;
    std::size_t worst_row = 0;
    std::size_t worst_column = 0;
    double worst_ad = 0.0;
    double worst_finite_difference = 0.0;

    const double step_factor = std::pow(
        std::numeric_limits<double>::epsilon(), 0.2);
    for (std::size_t fixture_index = 0;
         fixture_index < fixtures.size();
         ++fixture_index) {
        const DerivativeFixture& fixture =
            fixtures[fixture_index];
        check(
            "derivative fixture is in KS domain",
            fixture.metric.valid_point(fixture.point));

        const Vec3 analytic =
            fixture.metric.radial_coordinate_gradient(
                fixture.point);
        const Vec3 independent =
            independent_dual_radius_gradient(
                fixture.point,
                fixture.metric.spin_length());
        for (std::size_t component = 0;
             component < 3;
             ++component) {
            maximum_gradient_error = std::max(
                maximum_gradient_error,
                normalized_difference(
                    analytic[component],
                    independent[component]));
        }

        const auto ad =
            fixture.metric.contravariant_derivatives(
                fixture.point);
        Mat4 zero{};
        check(
            "stationary inverse derivative is exactly zero",
            ad[0] == zero);
        for (std::size_t coordinate = 1;
             coordinate < 4;
             ++coordinate) {
            const double scale = std::max(
                {fixture.metric.mass(),
                 std::fabs(fixture.point.v[coordinate]),
                 1.0});
            const Mat4 finite_difference =
                five_point_inverse_derivative(
                    fixture.metric,
                    fixture.point,
                    coordinate,
                    step_factor * scale);
            for (std::size_t row = 0; row < 4; ++row) {
                for (std::size_t column = 0;
                     column < 4;
                     ++column) {
                    const double error =
                        normalized_difference(
                            ad[coordinate][row][column],
                            finite_difference[row][column]);
                    if (error > maximum_derivative_error) {
                        maximum_derivative_error = error;
                        worst_fixture = fixture_index;
                        worst_coordinate = coordinate;
                        worst_row = row;
                        worst_column = column;
                        worst_ad =
                            ad[coordinate][row][column];
                        worst_finite_difference =
                            finite_difference[row][column];
                    }
                }
            }
        }
    }

    check(
        "analytic KS radius gradient matches independent Dual4",
        maximum_gradient_error < 1.0e-12,
        maximum_gradient_error);
    check(
        "KS AD inverse derivative matches five-point difference",
        maximum_derivative_error < 3.0e-8,
        maximum_derivative_error);
    std::cout << "  worst_derivative fixture=" << worst_fixture
              << " coordinate=" << worst_coordinate
              << " row=" << worst_row
              << " column=" << worst_column
              << " ad=" << worst_ad
              << " finite_difference="
              << worst_finite_difference << '\n';

    const Contravariant4 ring{
        Vec4{{0.0, positive.spin_length(), 0.0, 0.0}}};
    try {
        (void)positive.contravariant_derivatives(ring);
        check("invalid derivative point throws", false);
    } catch (const std::domain_error&) {
        check("invalid derivative point throws", true);
    }
    try {
        (void)positive.radial_coordinate_gradient(ring);
        check("invalid gradient point throws", false);
    } catch (const std::domain_error&) {
        check("invalid gradient point throws", true);
    }

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
