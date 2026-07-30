#include "solar/relativity/fluid_model.h"

#include "kerr_fluid_kinematics.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace solar::relativity {
namespace {

void require_positive_finite(
    double value,
    const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(message);
    }
}

} // namespace

AnalyticOpticallyThinTorus::
AnalyticOpticallyThinTorus(
    AnalyticOpticallyThinTorusConfig config)
    : config_(std::move(config)) {
    require_positive_finite(
        config_.mass_M,
        "torus mass must be finite and positive");
    if (!std::isfinite(config_.spin_chi) ||
        std::fabs(config_.spin_chi) >= 1.0) {
        throw std::invalid_argument(
            "torus spin must be finite with abs(chi)<1");
    }
    if (!is_valid_orbit_sense(config_.sense)) {
        throw std::invalid_argument(
            "torus orbit sense is not recognized");
    }
    require_positive_finite(
        config_.center_radius_M,
        "torus center radius must be finite and positive");
    require_positive_finite(
        config_.radial_width_M,
        "torus radial width must be finite and positive");
    require_positive_finite(
        config_.angular_width,
        "torus angular width must be finite and positive");
    if (config_.angular_width > 1.0) {
        throw std::invalid_argument(
            "torus angular width cannot exceed one");
    }
    require_positive_finite(
        config_.density_scale,
        "torus density scale must be finite and positive");
    require_positive_finite(
        config_.temperature_scale,
        "torus temperature scale must be finite and positive");
    if (!std::isfinite(config_.temperature_power) ||
        config_.temperature_power < 0.0) {
        throw std::invalid_argument(
            "torus temperature power must be finite and non-negative");
    }
    if (!std::isfinite(
            config_.density_cutoff_fraction) ||
        config_.density_cutoff_fraction <= 0.0 ||
        config_.density_cutoff_fraction >= 1.0) {
        throw std::invalid_argument(
            "torus cutoff fraction must lie strictly between zero and one");
    }
}

FluidSample AnalyticOpticallyThinTorus::sample(
    const Metric& metric,
    const Contravariant4& x) const {
    const detail::KerrFluidPoint point =
        detail::evaluate_kerr_circular_fluid_point(
            metric,
            x,
            config_.mass_M,
            config_.spin_chi,
            config_.sense);
    const double radial_offset =
        (point.radius - config_.center_radius_M) /
        config_.radial_width_M;
    const double angular_offset =
        point.equatorial_height /
        config_.angular_width;
    const double shape = std::exp(
        -0.5 *
        (radial_offset * radial_offset +
         angular_offset * angular_offset));
    if (!std::isfinite(shape)) {
        throw std::domain_error(
            "analytic torus shape is non-finite");
    }
    if (shape < config_.density_cutoff_fraction) {
        return {};
    }

    const double density =
        config_.density_scale * shape;
    const double temperature =
        config_.temperature_scale *
        std::pow(shape, config_.temperature_power);
    if (!std::isfinite(density) ||
        !std::isfinite(temperature)) {
        throw std::domain_error(
            "analytic torus profile is non-finite");
    }
    return FluidSample{
        true,
        density,
        temperature,
        point.caller_four_velocity,
    };
}

} // namespace solar::relativity
