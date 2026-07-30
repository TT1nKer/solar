#include "solar/relativity/observer.h"

#include "observer_validation.h"
#include "solar/relativity/spacetime_algebra.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace solar::relativity {
namespace {

ObserverResult observer_failure(
    ObserverError error,
    std::string message) {
    return ObserverResult{
        error, std::nullopt, std::move(message)};
}

double spatial_determinant(const Tetrad& tetrad) {
    const Vec4& first = tetrad.basis[1].v;
    const Vec4& second = tetrad.basis[2].v;
    const Vec4& third = tetrad.basis[3].v;
    return first[1] *
               (second[2] * third[3] -
                second[3] * third[2]) -
           first[2] *
               (second[1] * third[3] -
                second[3] * third[1]) +
           first[3] *
               (second[1] * third[2] -
                second[2] * third[1]);
}

Contravariant4 coordinate_seed(std::size_t coordinate) {
    Contravariant4 seed;
    seed.v[coordinate] = 1.0;
    return seed;
}

} // namespace

ObserverResult make_arbitrary_observer(
    const Metric& metric,
    const Contravariant4& x,
    const Contravariant4& four_velocity,
    const std::array<Contravariant4, 3>& spatial_seeds) {
    if (!x.v.all_finite() ||
        !four_velocity.v.all_finite()) {
        return observer_failure(
            ObserverError::NonFiniteInput,
            "observer position and four-velocity must be finite");
    }
    for (const Contravariant4& seed : spatial_seeds) {
        if (!seed.v.all_finite()) {
            return observer_failure(
                ObserverError::NonFiniteInput,
                "observer spatial seeds must be finite");
        }
    }
    if (!metric.valid_point(x)) {
        return observer_failure(
            ObserverError::InvalidMetricPoint,
            "observer position is outside the metric domain");
    }

    Mat4 covariant;
    try {
        covariant = metric.covariant(x);
    } catch (const std::domain_error&) {
        return observer_failure(
            ObserverError::InvalidMetricPoint,
            "observer metric evaluation failed");
    }

    const double velocity_norm = metric_inner_product(
        covariant, four_velocity, four_velocity);
    if (!std::isfinite(velocity_norm) ||
        std::fabs(velocity_norm + 1.0) >=
            detail::observer_tetrad_tolerance) {
        return observer_failure(
            ObserverError::FourVelocityNotUnitTimelike,
            "observer four-velocity must be unit timelike");
    }

    Tetrad tetrad;
    tetrad.basis[0] = four_velocity;
    for (std::size_t spatial = 0; spatial < 3; ++spatial) {
        Contravariant4 candidate = spatial_seeds[spatial];
        const double time_projection = metric_inner_product(
            covariant, candidate, tetrad.basis[0]);
        candidate.v = candidate.v +
                      time_projection * tetrad.basis[0].v;

        for (std::size_t completed = 0;
             completed < spatial;
             ++completed) {
            const Contravariant4& basis =
                tetrad.basis[completed + 1];
            const double projection = metric_inner_product(
                covariant, candidate, basis);
            candidate.v = candidate.v - projection * basis.v;
        }

        const double norm_squared = metric_inner_product(
            covariant, candidate, candidate);
        if (!std::isfinite(norm_squared) ||
            norm_squared <= 0.0) {
            return observer_failure(
                ObserverError::DegenerateSpatialSeed,
                "observer spatial seeds do not span the rest space");
        }
        candidate.v =
            candidate.v / std::sqrt(norm_squared);
        tetrad.basis[spatial + 1] = candidate;
    }

    const double determinant = spatial_determinant(tetrad);
    if (!std::isfinite(determinant) || determinant == 0.0) {
        return observer_failure(
            ObserverError::DegenerateSpatialSeed,
            "observer spatial basis has no coordinate orientation");
    }
    if (determinant < 0.0) {
        tetrad.basis[3].v = -tetrad.basis[3].v;
    }

    ObserverFrame observer{x, tetrad};
    if (tetrad_orthonormality_error(metric, observer) >=
        detail::observer_tetrad_tolerance) {
        return observer_failure(
            ObserverError::TetradValidationFailure,
            "observer tetrad exceeds the orthonormality tolerance");
    }
    return ObserverResult{
        ObserverError::None,
        std::move(observer),
        {},
    };
}

ObserverResult make_static_observer(
    const Metric& metric,
    const Contravariant4& x) {
    if (!x.v.all_finite()) {
        return observer_failure(
            ObserverError::NonFiniteInput,
            "static observer position must be finite");
    }
    if (!metric.valid_point(x)) {
        return observer_failure(
            ObserverError::InvalidMetricPoint,
            "static observer position is outside the metric domain");
    }

    Mat4 covariant;
    try {
        covariant = metric.covariant(x);
    } catch (const std::domain_error&) {
        return observer_failure(
            ObserverError::InvalidMetricPoint,
            "static observer metric evaluation failed");
    }
    if (!std::isfinite(covariant[0][0]) ||
        covariant[0][0] >= 0.0) {
        return observer_failure(
            ObserverError::StaticWorldlineNotTimelike,
            "static worldline is not timelike at this point");
    }

    Contravariant4 four_velocity;
    four_velocity.v[0] =
        1.0 / std::sqrt(-covariant[0][0]);
    const std::array<Contravariant4, 3> seeds{{
        coordinate_seed(1),
        coordinate_seed(2),
        coordinate_seed(3),
    }};
    return make_arbitrary_observer(
        metric, x, four_velocity, seeds);
}

ObserverResult make_look_at_observer(
    const Metric& metric,
    const Contravariant4& x,
    const Contravariant4& four_velocity,
    const LookAtAttitude& attitude) {
    if (!attitude.look_direction.v.all_finite() ||
        !attitude.up_reference.v.all_finite()) {
        return observer_failure(
            ObserverError::NonFiniteInput,
            "look-at attitude vectors must be finite");
    }

    for (std::size_t coordinate = 1;
         coordinate < 4;
         ++coordinate) {
        const std::array<Contravariant4, 3> seeds{{
            attitude.look_direction,
            attitude.up_reference,
            coordinate_seed(coordinate),
        }};
        ObserverResult result = make_arbitrary_observer(
            metric, x, four_velocity, seeds);
        if (result ||
            result.error !=
                ObserverError::DegenerateSpatialSeed) {
            return result;
        }
    }
    return observer_failure(
        ObserverError::DegenerateSpatialSeed,
        "look and up vectors do not define a spatial attitude");
}

} // namespace solar::relativity
