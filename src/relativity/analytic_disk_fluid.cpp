#include "solar/relativity/fluid_model.h"

#include "kerr_fluid_kinematics.h"
#include "solar/relativity/kerr_bl_metric.h"

#include <algorithm>
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

AnalyticCircularDiskFluid::AnalyticCircularDiskFluid(
    AnalyticCircularDiskConfig config)
    : config_(std::move(config)),
      inner_radius_M_(0.0) {
    require_positive_finite(
        config_.mass_M,
        "disk mass must be finite and positive");
    if (!std::isfinite(config_.spin_chi) ||
        std::fabs(config_.spin_chi) >= 1.0) {
        throw std::invalid_argument(
            "disk spin must be finite with abs(chi)<1");
    }
    if (!is_valid_orbit_sense(config_.sense)) {
        throw std::invalid_argument(
            "disk orbit sense is not recognized");
    }

    const KerrBoyerLindquistMetric metric(
        config_.mass_M, config_.spin_chi);
    inner_radius_M_ =
        config_.inner_radius_M.has_value()
            ? *config_.inner_radius_M
            : kerr_isco_radius(metric, config_.sense);
    require_positive_finite(
        inner_radius_M_,
        "disk inner radius must be finite and positive");
    if (inner_radius_M_ <=
        kerr_equatorial_photon_radius(
            metric, config_.sense)) {
        throw std::invalid_argument(
            "disk inner radius must exceed the circular photon orbit");
    }
    require_positive_finite(
        config_.outer_radius_M,
        "disk outer radius must be finite and positive");
    if (config_.outer_radius_M <= inner_radius_M_) {
        throw std::invalid_argument(
            "disk outer radius must exceed the inner radius");
    }
    require_positive_finite(
        config_.density_scale,
        "disk density scale must be finite and positive");
    require_positive_finite(
        config_.temperature_scale,
        "disk temperature scale must be finite and positive");
    if (!std::isfinite(config_.density_power) ||
        config_.density_power < 0.0) {
        throw std::invalid_argument(
            "disk density power must be finite and non-negative");
    }
    require_positive_finite(
        config_.surface_height_tolerance,
        "disk surface tolerance must be finite and positive");
    if (config_.surface_height_tolerance > 1.0) {
        throw std::invalid_argument(
            "disk surface tolerance cannot exceed one");
    }
}

FluidSample AnalyticCircularDiskFluid::sample(
    const Metric& metric,
    const Contravariant4& x) const {
    const detail::KerrFluidLocation location =
        detail::locate_kerr_fluid_point(
            metric,
            x,
            config_.mass_M,
            config_.spin_chi);
    if (location.radius < inner_radius_M_ ||
        location.radius > config_.outer_radius_M ||
        std::fabs(location.equatorial_height) >
            config_.surface_height_tolerance) {
        return {};
    }
    const Contravariant4 four_velocity =
        detail::evaluate_kerr_circular_four_velocity(
            metric,
            x,
            location,
            config_.mass_M,
            config_.spin_chi,
            config_.sense);

    const double normalized_radius =
        location.radius / inner_radius_M_;
    const double density =
        config_.density_scale *
        std::pow(
            normalized_radius,
            -config_.density_power);
    const double flux_shape =
        std::pow(normalized_radius, -3.0) *
        std::max(
            0.0,
            1.0 -
                std::sqrt(
                    1.0 / normalized_radius));
    const double temperature =
        config_.temperature_scale *
        std::pow(flux_shape, 0.25);
    if (!std::isfinite(density) ||
        !std::isfinite(temperature)) {
        throw std::domain_error(
            "analytic disk profile is non-finite");
    }

    return FluidSample{
        true,
        density,
        temperature,
        four_velocity,
    };
}

} // namespace solar::relativity
