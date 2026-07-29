#include "solar/relativity/dual4.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <utility>

using namespace solar::relativity;

namespace {

template <typename Scalar>
using ScalarMatrix4 = std::array<std::array<Scalar, 4>, 4>;

int passed = 0;
int failed = 0;

void check(const char* name, bool condition, double diagnostic = 0.0) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name
              << " (max normalized error=" << diagnostic << ")\n";
    condition ? ++passed : ++failed;
}

ScalarMatrix4<Dual4> independent_kerr_covariant_dual(
    const Contravariant4& point, double mass, double spin_chi) {
    const double spin_a = mass * spin_chi;
    const Dual4 radius = Dual4::variable(point.v[1], 1);
    const Dual4 theta = Dual4::variable(point.v[2], 2);
    const Dual4 sine_theta = sin(theta);
    const Dual4 cosine_theta = cos(theta);
    const Dual4 radius_squared = radius * radius;
    const double spin_squared = spin_a * spin_a;
    const Dual4 sigma =
        radius_squared + spin_squared * cosine_theta * cosine_theta;
    const Dual4 delta =
        radius_squared - 2.0 * mass * radius + spin_squared;
    const Dual4 radius_spin_sum = radius_squared + spin_squared;
    const Dual4 A =
        radius_spin_sum * radius_spin_sum -
        spin_squared * delta * sine_theta * sine_theta;
    const Dual4 sine_squared = sine_theta * sine_theta;

    ScalarMatrix4<Dual4> metric{};
    metric[0][0] = -(1.0 - 2.0 * mass * radius / sigma);
    metric[0][3] =
        -2.0 * mass * spin_a * radius * sine_squared / sigma;
    metric[3][0] = metric[0][3];
    metric[1][1] = sigma / delta;
    metric[2][2] = sigma;
    metric[3][3] = A * sine_squared / sigma;
    return metric;
}

double normalized_difference(double actual, double expected) {
    return std::fabs(actual - expected) /
           std::max({1.0, std::fabs(actual), std::fabs(expected)});
}

double kerr_dual_derivative_error(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& point) {
    const auto analytic = metric.covariant_derivatives(point);
    const auto dual = independent_kerr_covariant_dual(
        point, metric.mass(), metric.spin_chi());
    double maximum = 0.0;
    for (std::size_t coordinate = 0; coordinate < 4; ++coordinate) {
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0; column < 4; ++column) {
                maximum = std::max(
                    maximum,
                    normalized_difference(
                        analytic[coordinate][row][column],
                        dual[row][column].derivative[coordinate]));
            }
        }
    }
    return maximum;
}

Mat4 five_point_inverse_derivative(
    const Metric& metric,
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

double finite_difference_error(
    const Metric& metric,
    const Contravariant4& point) {
    const auto analytic = metric.contravariant_derivatives(point);
    double maximum = 0.0;
    for (std::size_t coordinate = 0; coordinate < 4; ++coordinate) {
        double step = 1.0e-5;
        if (coordinate == 1) {
            step = 2.0e-4 * std::max(1.0, std::fabs(point.v[1]));
        }
        const Mat4 finite_difference =
            five_point_inverse_derivative(
                metric, point, coordinate, step);
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t column = 0; column < 4; ++column) {
                maximum = std::max(
                    maximum,
                    normalized_difference(
                        analytic[coordinate][row][column],
                        finite_difference[row][column]));
            }
        }
    }
    return maximum;
}

} // namespace

int main() {
    double maximum_dual_error = 0.0;
    double maximum_kerr_fd_error = 0.0;

    for (const auto& test_case : {
             std::pair{
                 KerrBoyerLindquistMetric(1.0, 0.9),
                 Contravariant4{Vec4{{0.3, 10.0, 1.1, -0.2}}}},
             std::pair{
                 KerrBoyerLindquistMetric(2.3, -0.7),
                 Contravariant4{Vec4{{-4.0, 15.0, 2.0, 1.7}}}},
             std::pair{
                 KerrBoyerLindquistMetric(1.0, 0.2),
                 Contravariant4{Vec4{{2.0, 3.0, 0.8, 0.5}}}},
         }) {
        maximum_dual_error = std::max(
            maximum_dual_error,
            kerr_dual_derivative_error(
                test_case.first, test_case.second));
        maximum_kerr_fd_error = std::max(
            maximum_kerr_fd_error,
            finite_difference_error(
                test_case.first, test_case.second));
    }

    check("Kerr analytic covariant derivative matches Dual4",
          maximum_dual_error < 1.0e-12, maximum_dual_error);
    check("Kerr inverse derivative matches five-point difference",
          maximum_kerr_fd_error < 1.0e-8, maximum_kerr_fd_error);

    const SchwarzschildBoyerLindquistMetric schwarzschild(1.4);
    const Contravariant4 schwarzschild_point{
        Vec4{{1.0, 12.0, 1.2, -2.0}}};
    const double schwarzschild_fd_error =
        finite_difference_error(schwarzschild, schwarzschild_point);
    check("Schwarzschild inverse derivative matches five-point difference",
          schwarzschild_fd_error < 1.0e-8, schwarzschild_fd_error);

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
