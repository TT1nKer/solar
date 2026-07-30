#include "solar/relativity/kerr_orbits.h"

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>

namespace solar::relativity {
namespace {

ObserverResult circular_observer_failure(
    ObserverError error,
    std::string message) {
    return ObserverResult{
        error, std::nullopt, std::move(message)};
}

Contravariant4 coordinate_seed(std::size_t coordinate) {
    Contravariant4 seed;
    seed.v[coordinate] = 1.0;
    return seed;
}

} // namespace

ObserverResult make_equatorial_circular_observer(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& x,
    OrbitSense sense) {
    if (!is_valid_orbit_sense(sense)) {
        throw std::invalid_argument(
            "Kerr orbit sense is not recognized");
    }
    if (!x.v.all_finite()) {
        return circular_observer_failure(
            ObserverError::NonFiniteInput,
            "circular observer position must be finite");
    }
    if (std::fabs(
            x.v[2] - 1.5707963267948966) > 1.0e-12) {
        return circular_observer_failure(
            ObserverError::CircularWorldlineNotTimelike,
            "circular observer must lie in the equatorial plane");
    }
    if (!metric.valid_point(x)) {
        return circular_observer_failure(
            ObserverError::InvalidMetricPoint,
            "circular observer position is outside the Kerr BL domain");
    }

    const CircularOrbitResult orbit =
        evaluate_equatorial_circular_timelike_orbit(
            metric, x.v[1], sense);
    if (!orbit) {
        return circular_observer_failure(
            ObserverError::CircularWorldlineNotTimelike,
            orbit.message);
    }

    Mat4 covariant;
    try {
        covariant = metric.covariant(x);
    } catch (const std::domain_error&) {
        return circular_observer_failure(
            ObserverError::InvalidMetricPoint,
            "circular observer metric evaluation failed");
    }
    const double angular_velocity =
        orbit.orbit->angular_velocity;
    const double normalization_squared =
        -(covariant[0][0] +
          2.0 * angular_velocity * covariant[0][3] +
          angular_velocity * angular_velocity *
              covariant[3][3]);
    if (!std::isfinite(normalization_squared) ||
        normalization_squared <= 0.0) {
        return circular_observer_failure(
            ObserverError::CircularWorldlineNotTimelike,
            "circular worldline normalization is non-positive");
    }

    Contravariant4 four_velocity;
    four_velocity.v[0] =
        1.0 / std::sqrt(normalization_squared);
    four_velocity.v[3] =
        angular_velocity * four_velocity.v[0];
    const std::array<Contravariant4, 3> seeds{{
        coordinate_seed(1),
        coordinate_seed(2),
        coordinate_seed(3),
    }};
    return make_arbitrary_observer(
        metric, x, four_velocity, seeds);
}

} // namespace solar::relativity
