#include "solar/relativity/kerr_chart_transform.h"

#include "kerr_chart_fields.h"
#include "solar/constants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

constexpr double axis_sine_floor = 1.0e-12;
constexpr double jacobian_identity_tolerance = 1.0e-9;

void require_safe_boyer_lindquist(
    const KerrSchildCartesianMetric& metric,
    double minimum_radius,
    const Contravariant4& point) {
    if (!point.v.all_finite()) {
        throw std::domain_error(
            "Boyer-Lindquist transform point must be finite");
    }
    const double radius = point.v[1];
    const double theta = point.v[2];
    if (!(radius > minimum_radius) ||
        !(theta > 0.0 && theta < constants::PI) ||
        std::fabs(std::sin(theta)) <= axis_sine_floor) {
        throw std::domain_error(
            "Boyer-Lindquist point is outside the safe transform overlap");
    }
    const double delta =
        radius * radius -
        2.0 * metric.mass() * radius +
        metric.spin_length() * metric.spin_length();
    if (!std::isfinite(delta) || delta <= 0.0) {
        throw std::domain_error(
            "Boyer-Lindquist transform delta must be finite and positive");
    }
}

double matrix_identity_error(const Mat4& matrix) {
    double maximum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            maximum = std::max(
                maximum,
                std::fabs(matrix[row][column] - expected));
        }
    }
    return maximum;
}

Mat4 checked_inverse(const Mat4& jacobian) {
    if (!all_finite(jacobian)) {
        throw std::domain_error(
            "Kerr chart Jacobian is non-finite");
    }
    const Mat4 result = inverse(jacobian);
    const double identity_error = std::max(
        matrix_identity_error(multiply(result, jacobian)),
        matrix_identity_error(multiply(jacobian, result)));
    if (!std::isfinite(identity_error) ||
        identity_error > jacobian_identity_tolerance) {
        throw std::domain_error(
            "Kerr chart Jacobian is ill-conditioned");
    }
    return result;
}

double wrap_angle(double angle) {
    const double wrapped =
        std::remainder(angle, 2.0 * constants::PI);
    if (!std::isfinite(wrapped)) {
        throw std::domain_error(
            "Kerr chart azimuth is non-finite");
    }
    return wrapped;
}

std::array<double, 4> components(
    const Contravariant4& point) {
    return {
        point.v[0],
        point.v[1],
        point.v[2],
        point.v[3],
    };
}

} // namespace

KerrChartTransform::KerrChartTransform(
    double mass_M,
    double spin_chi,
    double overlap_margin_fraction)
    : kerr_schild_metric_(mass_M, spin_chi),
      overlap_margin_M_(
          mass_M * overlap_margin_fraction),
      minimum_overlap_radius_M_(
          kerr_schild_metric_.outer_horizon_radius() +
          overlap_margin_M_) {
    if (!std::isfinite(overlap_margin_fraction) ||
        overlap_margin_fraction <= 0.0 ||
        !std::isfinite(overlap_margin_M_) ||
        overlap_margin_M_ <= 0.0 ||
        !std::isfinite(minimum_overlap_radius_M_)) {
        throw std::invalid_argument(
            "Kerr chart overlap margin must be finite, positive, "
            "and representable");
    }
}

Contravariant4
KerrChartTransform::position_to_kerr_schild(
    const Contravariant4& boyer_lindquist) const {
    require_safe_boyer_lindquist(
        kerr_schild_metric_,
        minimum_overlap_radius_M_,
        boyer_lindquist);
    const auto transformed =
        detail::boyer_lindquist_to_kerr_schild_position(
            components(boyer_lindquist),
            kerr_schild_metric_.mass(),
            kerr_schild_metric_.spin_length(),
            kerr_schild_metric_.outer_horizon_radius(),
            kerr_schild_metric_.inner_horizon_radius());
    const Contravariant4 result{Vec4{{
        transformed[0],
        transformed[1],
        transformed[2],
        transformed[3],
    }}};
    if (!result.v.all_finite() ||
        !kerr_schild_metric_.valid_point(result)) {
        throw std::domain_error(
            "Boyer-Lindquist to Kerr-Schild position is invalid");
    }
    return result;
}

Contravariant4
KerrChartTransform::position_to_boyer_lindquist(
    const Contravariant4& kerr_schild) const {
    if (!kerr_schild.v.all_finite() ||
        !kerr_schild_metric_.valid_point(kerr_schild)) {
        throw std::domain_error(
            "Kerr-Schild transform point is outside the metric domain");
    }
    const double radius =
        kerr_schild_metric_.radial_coordinate(kerr_schild);
    if (!(radius > minimum_overlap_radius_M_)) {
        throw std::domain_error(
            "Kerr-Schild point is outside the safe transform overlap");
    }

    const double cosine_theta = kerr_schild.v[3] / radius;
    if (!(cosine_theta > -1.0 && cosine_theta < 1.0)) {
        throw std::domain_error(
            "Kerr-Schild point is on the unresolved polar axis");
    }
    const double theta = std::acos(cosine_theta);
    if (std::fabs(std::sin(theta)) <= axis_sine_floor) {
        throw std::domain_error(
            "Kerr-Schild point is on the unresolved polar axis");
    }

    const double spin_a =
        kerr_schild_metric_.spin_length();
    const double tilde_phi = std::atan2(
        radius * kerr_schild.v[2] -
            spin_a * kerr_schild.v[1],
        radius * kerr_schild.v[1] +
            spin_a * kerr_schild.v[2]);
    if (!std::isfinite(tilde_phi)) {
        throw std::domain_error(
            "Kerr-Schild azimuth recovery failed");
    }
    const double time_offset =
        detail::kerr_ingoing_time_offset(
            radius,
            kerr_schild_metric_.mass(),
            kerr_schild_metric_.outer_horizon_radius(),
            kerr_schild_metric_.inner_horizon_radius());
    const double azimuth_offset =
        detail::kerr_ingoing_azimuth_offset(
            radius,
            spin_a,
            kerr_schild_metric_.outer_horizon_radius(),
            kerr_schild_metric_.inner_horizon_radius());
    const Contravariant4 result{Vec4{{
        kerr_schild.v[0] - time_offset,
        radius,
        theta,
        wrap_angle(tilde_phi - azimuth_offset),
    }}};
    if (!result.v.all_finite()) {
        throw std::domain_error(
            "Kerr-Schild to Boyer-Lindquist position is invalid");
    }
    return result;
}

Mat4
KerrChartTransform::boyer_lindquist_to_kerr_schild_jacobian(
    const Contravariant4& boyer_lindquist) const {
    require_safe_boyer_lindquist(
        kerr_schild_metric_,
        minimum_overlap_radius_M_,
        boyer_lindquist);
    std::array<Dual4, 4> dual_coordinates{};
    for (std::size_t coordinate = 0;
         coordinate < 4;
         ++coordinate) {
        dual_coordinates[coordinate] =
            Dual4::variable(
                boyer_lindquist.v[coordinate],
                coordinate);
    }
    const auto transformed =
        detail::boyer_lindquist_to_kerr_schild_position(
            dual_coordinates,
            kerr_schild_metric_.mass(),
            kerr_schild_metric_.spin_length(),
            kerr_schild_metric_.outer_horizon_radius(),
            kerr_schild_metric_.inner_horizon_radius());
    Mat4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result[row][column] =
                transformed[row].derivative[column];
        }
    }
    (void)checked_inverse(result);
    return result;
}

Mat4
KerrChartTransform::kerr_schild_to_boyer_lindquist_jacobian(
    const Contravariant4& kerr_schild) const {
    const Contravariant4 boyer_lindquist =
        position_to_boyer_lindquist(kerr_schild);
    return checked_inverse(
        boyer_lindquist_to_kerr_schild_jacobian(
            boyer_lindquist));
}

} // namespace solar::relativity
