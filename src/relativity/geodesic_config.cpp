#include "geodesic_config_internal.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity {
namespace detail {

void validate_geodesic_config(
    const GeodesicIntegrationConfig& config) {
    if (!std::isfinite(config.initial_step) ||
        config.initial_step == 0.0) {
        throw std::invalid_argument(
            "geodesic initial step must be finite and non-zero");
    }
    if (!std::isfinite(config.min_step) ||
        config.min_step <= 0.0 ||
        !std::isfinite(config.max_step) ||
        config.max_step < config.min_step) {
        throw std::invalid_argument(
            "geodesic steps must satisfy 0 < min_step <= max_step");
    }
    if (config.max_rejections_per_step == 0 ||
        config.max_total_steps == 0) {
        throw std::invalid_argument(
            "geodesic rejection and total-step limits must be positive");
    }
    if (!std::isfinite(config.max_affine) ||
        config.max_affine <= 0.0) {
        throw std::invalid_argument(
            "geodesic affine limit must be positive and finite");
    }
    const auto valid_optional_limit = [](double limit) {
        return (std::isfinite(limit) && limit > 0.0) ||
               (std::isinf(limit) && limit > 0.0);
    };
    if (!valid_optional_limit(config.max_proper_time) ||
        !valid_optional_limit(config.max_coordinate_time)) {
        throw std::invalid_argument(
            "optional geodesic limits must be positive or infinity");
    }
    if (config.kind == GeodesicKind::Null &&
        std::isfinite(config.max_proper_time)) {
        throw std::invalid_argument(
            "null geodesics do not have a proper-time limit");
    }
    if (!std::isfinite(config.constraint_tolerance) ||
        config.constraint_tolerance <= 0.0) {
        throw std::invalid_argument(
            "geodesic constraint tolerance must be positive and finite");
    }
    if (!std::isfinite(config.dopri5.relative_tolerance) ||
        config.dopri5.relative_tolerance <= 0.0) {
        throw std::invalid_argument(
            "geodesic relative tolerance must be positive and finite");
    }
    for (const double tolerance :
         config.dopri5.absolute_tolerance) {
        if (!std::isfinite(tolerance) || tolerance <= 0.0) {
            throw std::invalid_argument(
                "geodesic absolute tolerances must be positive and finite");
        }
    }
    if (!std::isfinite(config.dopri5.safety) ||
        config.dopri5.safety <= 0.0 ||
        config.dopri5.safety >= 1.0 ||
        !std::isfinite(config.dopri5.min_factor) ||
        config.dopri5.min_factor <= 0.0 ||
        config.dopri5.min_factor > 1.0 ||
        !std::isfinite(config.dopri5.max_factor) ||
        config.dopri5.max_factor < 1.0 ||
        config.dopri5.min_factor >
            config.dopri5.max_factor) {
        throw std::invalid_argument(
            "geodesic DOPRI5 controller configuration is invalid");
    }
    if (config.dopri5.error_norm !=
            numerics::ErrorNorm::RootMeanSquare &&
        config.dopri5.error_norm !=
            numerics::ErrorNorm::Maximum) {
        throw std::invalid_argument(
            "geodesic DOPRI5 error norm is not recognized");
    }
}

} // namespace detail

GeodesicIntegrationConfig
GeodesicIntegrationConfig::cpu_reference(
    GeodesicKind kind,
    double mass_scale,
    double initial_step,
    double max_step,
    double max_affine) {
    if (!std::isfinite(mass_scale) || mass_scale <= 0.0) {
        throw std::invalid_argument(
            "geodesic mass scale must be positive and finite");
    }

    numerics::Dopri5Config<8> dopri5{};
    dopri5.absolute_tolerance[0] = 1.0e-11 * mass_scale;
    for (std::size_t component = 1; component < 4; ++component) {
        dopri5.absolute_tolerance[component] =
            1.0e-12 * mass_scale;
    }
    for (std::size_t component = 4; component < 8; ++component) {
        dopri5.absolute_tolerance[component] = 1.0e-12;
    }
    dopri5.relative_tolerance = 1.0e-11;
    dopri5.safety = 0.9;
    dopri5.min_factor = 0.2;
    dopri5.max_factor = 5.0;
    dopri5.error_norm = numerics::ErrorNorm::RootMeanSquare;

    GeodesicIntegrationConfig config{
        kind,
        dopri5,
        initial_step,
        1.0e-12 * mass_scale,
        max_step,
        12,
        2000000,
        max_affine,
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        1.0e-10,
        false,
        false,
    };
    detail::validate_geodesic_config(config);
    return config;
}

} // namespace solar::relativity
