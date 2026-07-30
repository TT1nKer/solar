#include "solar/relativity/local_initialization.h"

#include "observer_validation.h"
#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace solar::relativity {
namespace {

constexpr double constraint_tolerance = 1.0e-10;
constexpr double frequency_tolerance = 1.0e-12;

InitialStateResult initialization_failure(
    InitialStateError error,
    std::string message) {
    return InitialStateResult{
        error,
        std::nullopt,
        std::numeric_limits<double>::quiet_NaN(),
        std::move(message),
    };
}

InitialStateResult finish_initialization(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec4& local_vector,
    double affine,
    GeodesicKind kind,
    bool normalize_frequency) {
    if (!std::isfinite(affine) ||
        !local_vector.all_finite()) {
        return initialization_failure(
            InitialStateError::NonFiniteInput,
            "local state and affine parameter must be finite");
    }
    if (!metric.valid_point(observer.x) ||
        tetrad_orthonormality_error(metric, observer) >=
            detail::observer_tetrad_tolerance) {
        return initialization_failure(
            InitialStateError::InvalidObserverFrame,
            "observer frame is invalid or non-orthonormal");
    }

    Mat4 covariant;
    try {
        covariant = metric.covariant(observer.x);
    } catch (const std::domain_error&) {
        return initialization_failure(
            InitialStateError::InvalidObserverFrame,
            "observer metric evaluation failed");
    }

    const Contravariant4 coordinate_vector =
        tetrad_to_coordinate(observer.tetrad, local_vector);
    if (!coordinate_vector.v.all_finite()) {
        return initialization_failure(
            InitialStateError::NonFiniteInput,
            "tetrad expansion produced a non-finite vector");
    }

    Covariant4 momentum =
        lower_index(covariant, coordinate_vector);
    double frequency = observer_measured_frequency(
        momentum, observer.tetrad.basis[0]);
    if (!std::isfinite(frequency) || frequency <= 0.0) {
        return initialization_failure(
            InitialStateError::NonFutureDirected,
            "local state is not future-directed for this observer");
    }

    if (normalize_frequency) {
        momentum.v = momentum.v / frequency;
        frequency = observer_measured_frequency(
            momentum, observer.tetrad.basis[0]);
        if (!std::isfinite(frequency) ||
            std::fabs(frequency - 1.0) >
                frequency_tolerance) {
            return initialization_failure(
                InitialStateError::NonFutureDirected,
                "photon frequency normalization failed");
        }
    }
    if (!momentum.v.all_finite()) {
        return initialization_failure(
            InitialStateError::NonFiniteInput,
            "canonical momentum is non-finite");
    }

    const PhaseSpaceState state{
        affine, observer.x, momentum};
    double constraint;
    try {
        constraint = hamiltonian_constraint_error(
            metric, state, kind);
    } catch (const std::domain_error&) {
        return initialization_failure(
            InitialStateError::ConstraintViolation,
            "Hamiltonian constraint evaluation failed");
    }
    if (!std::isfinite(constraint) ||
        constraint > constraint_tolerance) {
        return initialization_failure(
            InitialStateError::ConstraintViolation,
            "initialized state exceeds the Hamiltonian constraint gate");
    }

    return InitialStateResult{
        InitialStateError::None,
        state,
        frequency,
        {},
    };
}

} // namespace

double observer_measured_frequency(
    const Covariant4& momentum,
    const Contravariant4& observer_velocity) noexcept {
    return -covector_vector_pairing(
        momentum, observer_velocity);
}

InitialStateResult initialize_local_photon(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec3& local_direction,
    double affine) {
    if (!local_direction.all_finite() ||
        !std::isfinite(affine)) {
        return initialization_failure(
            InitialStateError::NonFiniteInput,
            "photon direction and affine parameter must be finite");
    }

    const double direction_norm = std::hypot(
        local_direction[0],
        local_direction[1],
        local_direction[2]);
    if (!std::isfinite(direction_norm) ||
        direction_norm == 0.0) {
        return initialization_failure(
            InitialStateError::InvalidLocalDirection,
            "photon direction must be finite and nonzero");
    }

    const Vec4 local_photon{{
        1.0,
        local_direction[0] / direction_norm,
        local_direction[1] / direction_norm,
        local_direction[2] / direction_norm,
    }};
    return finish_initialization(
        metric,
        observer,
        local_photon,
        affine,
        GeodesicKind::Null,
        true);
}

InitialStateResult initialize_local_timelike(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec3& local_velocity,
    double affine) {
    if (!local_velocity.all_finite() ||
        !std::isfinite(affine)) {
        return initialization_failure(
            InitialStateError::NonFiniteInput,
            "local velocity and affine parameter must be finite");
    }

    const double speed = std::hypot(
        local_velocity[0],
        local_velocity[1],
        local_velocity[2]);
    if (!std::isfinite(speed) || speed >= 1.0) {
        return initialization_failure(
            InitialStateError::SuperluminalLocalVelocity,
            "local timelike speed must satisfy magnitude below one");
    }

    const double gamma =
        1.0 / std::sqrt(1.0 - speed * speed);
    const Vec4 local_timelike{{
        gamma,
        gamma * local_velocity[0],
        gamma * local_velocity[1],
        gamma * local_velocity[2],
    }};
    return finish_initialization(
        metric,
        observer,
        local_timelike,
        affine,
        GeodesicKind::TimelikeUnitMass,
        false);
}

} // namespace solar::relativity
