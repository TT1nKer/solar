#include "solar/relativity/observer.h"

#include "observer_validation.h"
#include "solar/relativity/kerr_bl_metric.h"

#include <cmath>
#include <optional>
#include <string>
#include <utility>

namespace solar::relativity {
namespace {

ObserverResult zamo_failure(
    ObserverError error,
    std::string message) {
    return ObserverResult{
        error, std::nullopt, std::move(message)};
}

} // namespace

ObserverResult make_zamo_observer(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& x) {
    if (!x.v.all_finite()) {
        return zamo_failure(
            ObserverError::NonFiniteInput,
            "ZAMO position must be finite");
    }
    if (!metric.valid_point(x)) {
        return zamo_failure(
            ObserverError::InvalidMetricPoint,
            "ZAMO position is outside the Kerr BL domain");
    }

    const double mass = metric.mass();
    const double spin = metric.spin_length();
    const double radius = x.v[1];
    const double theta = x.v[2];
    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    const double radius_squared = radius * radius;
    const double spin_squared = spin * spin;
    const double sigma =
        radius_squared +
        spin_squared * cos_theta * cos_theta;
    const double delta =
        radius_squared - 2.0 * mass * radius +
        spin_squared;
    const double radius_spin = radius_squared + spin_squared;
    const double area =
        radius_spin * radius_spin -
        spin_squared * delta * sin_theta * sin_theta;

    if (!std::isfinite(sigma) || sigma <= 0.0 ||
        !std::isfinite(delta) || delta <= 0.0 ||
        !std::isfinite(area) || area <= 0.0 ||
        !std::isfinite(sin_theta) || sin_theta == 0.0) {
        return zamo_failure(
            ObserverError::InvalidMetricPoint,
            "ZAMO factors are outside their finite exterior domain");
    }

    const double lapse_squared = sigma * delta / area;
    if (!std::isfinite(lapse_squared) ||
        lapse_squared <= 0.0) {
        return zamo_failure(
            ObserverError::InvalidMetricPoint,
            "ZAMO lapse is not positive and finite");
    }
    const double lapse = std::sqrt(lapse_squared);
    const double frame_dragging =
        2.0 * mass * spin * radius / area;

    Tetrad tetrad;
    tetrad.basis[0].v[0] = 1.0 / lapse;
    tetrad.basis[0].v[3] = frame_dragging / lapse;
    tetrad.basis[1].v[1] = std::sqrt(delta / sigma);
    tetrad.basis[2].v[2] = 1.0 / std::sqrt(sigma);
    tetrad.basis[3].v[3] =
        std::sqrt(sigma / area) / sin_theta;

    ObserverFrame observer{x, tetrad};
    if (!tetrad.basis[0].v.all_finite() ||
        !tetrad.basis[1].v.all_finite() ||
        !tetrad.basis[2].v.all_finite() ||
        !tetrad.basis[3].v.all_finite() ||
        tetrad_orthonormality_error(metric, observer) >=
            detail::observer_tetrad_tolerance) {
        return zamo_failure(
            ObserverError::TetradValidationFailure,
            "ZAMO tetrad exceeds the orthonormality tolerance");
    }

    return ObserverResult{
        ObserverError::None,
        std::move(observer),
        {},
    };
}

} // namespace solar::relativity
