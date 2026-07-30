#include "solar/relativity/kerr_bl_metric.h"

#include "solar/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity {
namespace {

constexpr double axis_sine_floor = 1.0e-12;

struct KerrScalars {
    double radius;
    double sine_theta;
    double cosine_theta;
    double sigma;
    double delta;
    double A;
    double dr_sigma;
    double dtheta_sigma;
    double dr_delta;
    double dr_A;
    double dtheta_A;
};

KerrScalars evaluate_scalars(
    const Contravariant4& x, double mass, double spin_a) {
    const double radius = x.v[1];
    const double sine_theta = std::sin(x.v[2]);
    const double cosine_theta = std::cos(x.v[2]);
    const double radius_squared = radius * radius;
    const double spin_squared = spin_a * spin_a;
    const double sine_squared = sine_theta * sine_theta;
    const double sigma =
        radius_squared + spin_squared * cosine_theta * cosine_theta;
    const double delta =
        radius_squared - 2.0 * mass * radius + spin_squared;
    const double radius_spin_sum = radius_squared + spin_squared;
    const double A =
        radius_spin_sum * radius_spin_sum -
        spin_squared * delta * sine_squared;
    const double dr_sigma = 2.0 * radius;
    const double dtheta_sigma =
        -2.0 * spin_squared * sine_theta * cosine_theta;
    const double dr_delta = 2.0 * radius - 2.0 * mass;
    const double dr_A =
        4.0 * radius * radius_spin_sum -
        spin_squared * dr_delta * sine_squared;
    const double dtheta_A =
        -2.0 * spin_squared * delta * sine_theta * cosine_theta;
    return {
        radius,
        sine_theta,
        cosine_theta,
        sigma,
        delta,
        A,
        dr_sigma,
        dtheta_sigma,
        dr_delta,
        dr_A,
        dtheta_A,
    };
}

void require_valid(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& x) {
    if (!metric.valid_point(x)) {
        throw std::domain_error(
            "point is outside the regular exterior Kerr BL domain");
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

KerrBoyerLindquistMetric::KerrBoyerLindquistMetric(
    double mass_M,
    double spin_chi,
    double horizon_margin_fraction)
    : mass_M_(mass_M),
      spin_chi_(spin_chi),
      spin_a_M_(mass_M * spin_chi),
      horizon_margin_M_(mass_M * horizon_margin_fraction),
      outer_horizon_M_(0.0),
      inner_horizon_M_(0.0) {
    if (!std::isfinite(mass_M) || mass_M <= 0.0) {
        throw std::invalid_argument(
            "Kerr mass must be finite and positive");
    }
    if (!std::isfinite(spin_chi) || std::fabs(spin_chi) >= 1.0) {
        throw std::invalid_argument(
            "Kerr dimensionless spin must be finite with abs(chi)<1");
    }
    if (!std::isfinite(horizon_margin_fraction) ||
        horizon_margin_fraction <= 0.0) {
        throw std::invalid_argument(
            "horizon margin fraction must be finite and positive");
    }
    if (!std::isfinite(spin_a_M_) ||
        !std::isfinite(horizon_margin_M_) ||
        horizon_margin_M_ <= 0.0) {
        throw std::invalid_argument(
            "Kerr derived length scales must remain finite and positive");
    }

    const double spin_squared = spin_chi_ * spin_chi_;
    const double horizon_factor = std::sqrt(1.0 - spin_squared);
    outer_horizon_M_ = mass_M_ * (1.0 + horizon_factor);
    inner_horizon_M_ =
        mass_M_ * spin_squared / (1.0 + horizon_factor);
    if (!std::isfinite(outer_horizon_M_) ||
        !std::isfinite(inner_horizon_M_) ||
        !std::isfinite(outer_horizon_M_ + horizon_margin_M_)) {
        throw std::invalid_argument(
            "Kerr horizon radii must remain finite");
    }
}

Mat4 KerrBoyerLindquistMetric::covariant(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const KerrScalars q = evaluate_scalars(x, mass_M_, spin_a_M_);
    const double sine_squared = q.sine_theta * q.sine_theta;

    Mat4 result{};
    result[0][0] =
        -(1.0 - 2.0 * mass_M_ * q.radius / q.sigma);
    result[0][3] =
        -2.0 * mass_M_ * spin_a_M_ * q.radius *
        sine_squared / q.sigma;
    result[3][0] = result[0][3];
    result[1][1] = q.sigma / q.delta;
    result[2][2] = q.sigma;
    result[3][3] = q.A * sine_squared / q.sigma;
    return result;
}

Mat4 KerrBoyerLindquistMetric::contravariant(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const KerrScalars q = evaluate_scalars(x, mass_M_, spin_a_M_);
    const double sine_squared = q.sine_theta * q.sine_theta;
    const double sigma_delta = q.sigma * q.delta;

    Mat4 result{};
    result[0][0] = -q.A / sigma_delta;
    result[0][3] =
        -2.0 * mass_M_ * spin_a_M_ * q.radius / sigma_delta;
    result[3][0] = result[0][3];
    result[1][1] = q.delta / q.sigma;
    result[2][2] = 1.0 / q.sigma;
    result[3][3] =
        (q.delta - spin_a_M_ * spin_a_M_ * sine_squared) /
        (sigma_delta * sine_squared);
    return result;
}

std::array<Mat4, 4>
KerrBoyerLindquistMetric::covariant_derivatives(
    const Contravariant4& x) const {
    require_valid(*this, x);
    const KerrScalars q = evaluate_scalars(x, mass_M_, spin_a_M_);
    const double sine_squared = q.sine_theta * q.sine_theta;
    const double sigma_squared = q.sigma * q.sigma;
    const double delta_squared = q.delta * q.delta;

    std::array<Mat4, 4> result{};
    Mat4& dr = result[1];
    Mat4& dtheta = result[2];

    dr[0][0] =
        2.0 * mass_M_ / q.sigma -
        2.0 * mass_M_ * q.radius * q.dr_sigma / sigma_squared;
    dtheta[0][0] =
        -2.0 * mass_M_ * q.radius * q.dtheta_sigma /
        sigma_squared;

    dr[0][3] =
        -2.0 * mass_M_ * spin_a_M_ * sine_squared / q.sigma +
        2.0 * mass_M_ * spin_a_M_ * q.radius *
            sine_squared * q.dr_sigma / sigma_squared;
    dr[3][0] = dr[0][3];
    dtheta[0][3] =
        -4.0 * mass_M_ * spin_a_M_ * q.radius *
            q.sine_theta * q.cosine_theta / q.sigma +
        2.0 * mass_M_ * spin_a_M_ * q.radius *
            sine_squared * q.dtheta_sigma / sigma_squared;
    dtheta[3][0] = dtheta[0][3];

    dr[1][1] =
        (q.dr_sigma * q.delta - q.sigma * q.dr_delta) /
        delta_squared;
    dtheta[1][1] = q.dtheta_sigma / q.delta;

    dr[2][2] = q.dr_sigma;
    dtheta[2][2] = q.dtheta_sigma;

    dr[3][3] =
        q.dr_A * sine_squared / q.sigma -
        q.A * sine_squared * q.dr_sigma / sigma_squared;
    dtheta[3][3] =
        (q.dtheta_A * sine_squared +
         2.0 * q.A * q.sine_theta * q.cosine_theta) /
            q.sigma -
        q.A * sine_squared * q.dtheta_sigma / sigma_squared;
    return result;
}

std::array<Mat4, 4>
KerrBoyerLindquistMetric::contravariant_derivatives(
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

bool KerrBoyerLindquistMetric::valid_point(
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
    if (!(radius > outer_horizon_M_ + horizon_margin_M_)) {
        return false;
    }

    const KerrScalars q = evaluate_scalars(x, mass_M_, spin_a_M_);
    if (!std::isfinite(q.sigma) || q.sigma <= 0.0 ||
        !std::isfinite(q.delta) ||
        !std::isfinite(q.A) || q.A <= 0.0) {
        return false;
    }
    const double scale = std::max(
        {std::fabs(radius), mass_M_, std::fabs(spin_a_M_)});
    const double delta_floor =
        64.0 * std::numeric_limits<double>::epsilon() * scale * scale;
    return q.delta > delta_floor;
}

double KerrBoyerLindquistMetric::outer_stationary_limit_radius(
    double theta) const {
    if (!std::isfinite(theta)) {
        throw std::invalid_argument(
            "stationary-limit theta must be finite");
    }
    const double cosine_theta = std::cos(theta);
    return mass_M_ * (
        1.0 + std::sqrt(
            1.0 -
            spin_chi_ * spin_chi_ *
                cosine_theta * cosine_theta));
}

} // namespace solar::relativity
