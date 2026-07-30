#include "kerr_separated_config_internal.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity {
namespace detail {

void validate_kerr_separated_config(
    const KerrSeparatedConfig& config) {
    if (!is_valid_geodesic_kind(config.kind)) {
        throw std::invalid_argument(
            "Kerr separated geodesic kind is not recognized");
    }
    if (!std::isfinite(config.initial_mino_step) ||
        config.initial_mino_step == 0.0) {
        throw std::invalid_argument(
            "initial Mino step must be finite and non-zero");
    }
    if (!std::isfinite(config.min_mino_step) ||
        config.min_mino_step <= 0.0 ||
        !std::isfinite(config.max_mino_step) ||
        config.max_mino_step < config.min_mino_step) {
        throw std::invalid_argument(
            "Mino steps must satisfy 0 < min <= max");
    }
    if (config.max_rejections_per_step == 0 ||
        config.max_total_steps == 0) {
        throw std::invalid_argument(
            "Kerr separated step limits must be positive");
    }
    if (!std::isfinite(config.max_affine) ||
        config.max_affine <= 0.0) {
        throw std::invalid_argument(
            "maximum affine displacement must be positive and finite");
    }
    if (!((std::isfinite(config.max_coordinate_time) &&
           config.max_coordinate_time > 0.0) ||
          (std::isinf(config.max_coordinate_time) &&
           config.max_coordinate_time > 0.0))) {
        throw std::invalid_argument(
            "coordinate-time limit must be positive or infinity");
    }
    if (!std::isfinite(config.potential_tolerance) ||
        config.potential_tolerance <= 0.0 ||
        !std::isfinite(config.root_tolerance) ||
        config.root_tolerance <= 0.0 ||
        !std::isfinite(config.critical_derivative_tolerance) ||
        config.critical_derivative_tolerance <= 0.0 ||
        !std::isfinite(config.polar_axis_tolerance) ||
        config.polar_axis_tolerance <= 0.0 ||
        config.polar_axis_tolerance >= 1.0) {
        throw std::invalid_argument(
            "Kerr separated tolerances are invalid");
    }
    if (!std::isfinite(config.dopri5.relative_tolerance) ||
        config.dopri5.relative_tolerance <= 0.0) {
        throw std::invalid_argument(
            "Kerr separated relative tolerance is invalid");
    }
    for (const double tolerance :
         config.dopri5.absolute_tolerance) {
        if (!std::isfinite(tolerance) ||
            tolerance <= 0.0) {
            throw std::invalid_argument(
                "Kerr separated absolute tolerance is invalid");
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
            "Kerr separated DOPRI5 controller is invalid");
    }
    if (config.dopri5.error_norm !=
            numerics::ErrorNorm::RootMeanSquare &&
        config.dopri5.error_norm !=
            numerics::ErrorNorm::Maximum) {
        throw std::invalid_argument(
            "Kerr separated DOPRI5 norm is not recognized");
    }
}

} // namespace detail

KerrSeparatedConfig KerrSeparatedConfig::cpu_reference(
    GeodesicKind kind,
    double mass_scale,
    double initial_mino_step,
    double max_mino_step,
    double max_affine) {
    if (!std::isfinite(mass_scale) ||
        mass_scale <= 0.0) {
        throw std::invalid_argument(
            "Kerr separated mass scale must be positive and finite");
    }

    numerics::Dopri5Config<5> dopri5{};
    dopri5.absolute_tolerance = {{
        1.0e-12 * mass_scale,
        1.0e-12 * mass_scale,
        1.0e-13,
        1.0e-13,
        1.0e-12 * mass_scale,
    }};
    dopri5.relative_tolerance = 1.0e-11;
    dopri5.safety = 0.9;
    dopri5.min_factor = 0.2;
    dopri5.max_factor = 5.0;
    dopri5.error_norm =
        numerics::ErrorNorm::RootMeanSquare;

    KerrSeparatedConfig config{
        kind,
        dopri5,
        initial_mino_step,
        1.0e-14 / mass_scale,
        max_mino_step,
        32,
        1'000'000,
        max_affine,
        std::numeric_limits<double>::infinity(),
        1.0e-12,
        1.0e-12,
        1.0e-10,
        1.0e-12,
    };
    detail::validate_kerr_separated_config(config);
    return config;
}

} // namespace solar::relativity
