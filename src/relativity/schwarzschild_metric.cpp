#include "solar/relativity/schwarzschild_metric.h"

#include "solar/constants.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

constexpr double axis_sine_floor = 1.0e-12;

void require_valid(
    const SchwarzschildBoyerLindquistMetric& metric,
    const Contravariant4& x) {
    if (!metric.valid_point(x)) {
        throw std::domain_error(
            "point is outside the regular exterior Schwarzschild BL domain");
    }
}

Mat4 negated(const Mat4& matrix) {
    Mat4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result[row][column] = -matrix[row][column];
        }
    }
    return result;
}

} // namespace

SchwarzschildBoyerLindquistMetric::
SchwarzschildBoyerLindquistMetric(
    double mass_M, double horizon_margin_fraction)
    : mass_M_(mass_M),
      horizon_margin_M_(horizon_margin_fraction * mass_M) {
    if (!std::isfinite(mass_M) || mass_M <= 0.0) {
        throw std::invalid_argument(
            "Schwarzschild mass must be finite and positive");
    }
    if (!std::isfinite(horizon_margin_fraction) ||
        horizon_margin_fraction <= 0.0) {
        throw std::invalid_argument(
            "horizon margin fraction must be finite and positive");
    }
}

Mat4 SchwarzschildBoyerLindquistMetric::covariant(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const double radius = x.v[1];
    const double sine_theta = std::sin(x.v[2]);
    const double f = 1.0 - 2.0 * mass_M_ / radius;
    const double radius_squared = radius * radius;

    Mat4 result{};
    result[0][0] = -f;
    result[1][1] = 1.0 / f;
    result[2][2] = radius_squared;
    result[3][3] = radius_squared * sine_theta * sine_theta;
    return result;
}

Mat4 SchwarzschildBoyerLindquistMetric::contravariant(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const double radius = x.v[1];
    const double sine_theta = std::sin(x.v[2]);
    const double f = 1.0 - 2.0 * mass_M_ / radius;
    const double radius_squared = radius * radius;

    Mat4 result{};
    result[0][0] = -1.0 / f;
    result[1][1] = f;
    result[2][2] = 1.0 / radius_squared;
    result[3][3] =
        1.0 / (radius_squared * sine_theta * sine_theta);
    return result;
}

std::array<Mat4, 4>
SchwarzschildBoyerLindquistMetric::covariant_derivatives(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const double radius = x.v[1];
    const double theta = x.v[2];
    const double sine_theta = std::sin(theta);
    const double cosine_theta = std::cos(theta);
    const double radius_squared = radius * radius;
    const double f = 1.0 - 2.0 * mass_M_ / radius;

    std::array<Mat4, 4> result{};
    result[1][0][0] = -2.0 * mass_M_ / radius_squared;
    result[1][1][1] =
        -2.0 * mass_M_ / (radius_squared * f * f);
    result[1][2][2] = 2.0 * radius;
    result[1][3][3] = 2.0 * radius * sine_theta * sine_theta;
    result[2][3][3] =
        2.0 * radius_squared * sine_theta * cosine_theta;
    return result;
}

std::array<Mat4, 4>
SchwarzschildBoyerLindquistMetric::contravariant_derivatives(
    const Contravariant4& x) const {
    const Mat4 inverse_metric = contravariant(x);
    const auto covariant_partials = covariant_derivatives(x);
    std::array<Mat4, 4> result{};
    for (std::size_t coordinate = 0; coordinate < 4; ++coordinate) {
        result[coordinate] = negated(multiply(
            multiply(inverse_metric, covariant_partials[coordinate]),
            inverse_metric));
    }
    return result;
}

bool SchwarzschildBoyerLindquistMetric::valid_point(
    const Contravariant4& x) const noexcept {
    if (!x.v.all_finite()) {
        return false;
    }
    const double radius = x.v[1];
    const double theta = x.v[2];
    if (!(theta > 0.0 && theta < constants::PI)) {
        return false;
    }
    if (std::fabs(std::sin(theta)) <= axis_sine_floor) {
        return false;
    }
    return radius > outer_horizon_radius() + horizon_margin_M_;
}

} // namespace solar::relativity
