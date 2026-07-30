#include "kerr_fluid_kinematics.h"

#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

double rotation_sign(
    double spin_chi,
    OrbitSense sense) {
    const double spin_sign =
        spin_chi < 0.0 ? -1.0 : 1.0;
    return sense == OrbitSense::Prograde
               ? spin_sign
               : -spin_sign;
}

void require_matching_parameters(
    double actual_mass,
    double actual_spin,
    double expected_mass,
    double expected_spin) {
    if (actual_mass != expected_mass ||
        actual_spin != expected_spin) {
        throw std::domain_error(
            "fluid model and Kerr metric parameters do not match");
    }
}

Contravariant4 apply_jacobian(
    const Mat4& jacobian,
    const Contravariant4& vector) {
    Contravariant4 result{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            result.v[row] +=
                jacobian[row][column] * vector.v[column];
        }
    }
    return result;
}

} // namespace

KerrFluidPoint evaluate_kerr_circular_fluid_point(
    const Metric& metric,
    const Contravariant4& x,
    double mass_M,
    double spin_chi,
    OrbitSense sense) {
    if (!is_valid_orbit_sense(sense)) {
        throw std::invalid_argument(
            "Kerr fluid orbit sense is not recognized");
    }
    if (!x.v.all_finite() || !metric.valid_point(x)) {
        throw std::domain_error(
            "Kerr fluid point is outside the metric domain");
    }

    const auto* bl_metric =
        dynamic_cast<const KerrBoyerLindquistMetric*>(&metric);
    const auto* ks_metric =
        dynamic_cast<const KerrSchildCartesianMetric*>(&metric);
    if (bl_metric == nullptr && ks_metric == nullptr) {
        throw std::domain_error(
            "analytic Kerr fluid requires a Kerr BL or KS metric");
    }

    KerrBoyerLindquistMetric model_bl_metric(
        mass_M, spin_chi);
    Contravariant4 bl_position;
    bool transform_to_ks = false;
    if (bl_metric != nullptr) {
        require_matching_parameters(
            bl_metric->mass(),
            bl_metric->spin_chi(),
            mass_M,
            spin_chi);
        bl_position = x;
    } else {
        require_matching_parameters(
            ks_metric->mass(),
            ks_metric->spin_chi(),
            mass_M,
            spin_chi);
        const KerrChartTransform transform(
            mass_M, spin_chi);
        bl_position =
            transform.position_to_boyer_lindquist(x);
        transform_to_ks = true;
    }

    if (!model_bl_metric.valid_point(bl_position)) {
        throw std::domain_error(
            "Kerr fluid BL point is outside the metric domain");
    }
    const double radius = bl_position.v[1];
    const double normalized_radius = radius / mass_M;
    const double sign = rotation_sign(spin_chi, sense);
    const double angular_velocity =
        sign /
        (mass_M *
         (normalized_radius *
              std::sqrt(normalized_radius) +
          sign * spin_chi));
    const Mat4 bl_covariant =
        model_bl_metric.covariant(bl_position);
    const double normalization_squared =
        -(bl_covariant[0][0] +
          2.0 * angular_velocity *
              bl_covariant[0][3] +
          angular_velocity * angular_velocity *
              bl_covariant[3][3]);
    if (!std::isfinite(angular_velocity) ||
        !std::isfinite(normalization_squared) ||
        normalization_squared <= 0.0) {
        throw std::domain_error(
            "analytic Kerr circular flow is not timelike");
    }

    Contravariant4 bl_velocity{};
    bl_velocity.v[0] =
        1.0 / std::sqrt(normalization_squared);
    bl_velocity.v[3] =
        angular_velocity * bl_velocity.v[0];
    Contravariant4 caller_velocity = bl_velocity;
    if (transform_to_ks) {
        const KerrChartTransform transform(
            mass_M, spin_chi);
        caller_velocity = apply_jacobian(
            transform
                .boyer_lindquist_to_kerr_schild_jacobian(
                    bl_position),
            bl_velocity);
    }

    const double caller_norm = metric_inner_product(
        metric.covariant(x),
        caller_velocity,
        caller_velocity);
    if (!caller_velocity.v.all_finite() ||
        !std::isfinite(caller_norm) ||
        std::fabs(caller_norm + 1.0) > 1.0e-10 ||
        caller_velocity.v[0] <= 0.0) {
        throw std::domain_error(
            "analytic Kerr circular velocity failed validation");
    }

    return KerrFluidPoint{
        bl_position,
        caller_velocity,
        radius,
        bl_position.v[2],
        std::cos(bl_position.v[2]),
    };
}

} // namespace solar::relativity::detail
