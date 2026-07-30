#pragma once

#include "solar/relativity/observer.h"

#include <limits>
#include <optional>
#include <string>

namespace solar::relativity {

enum class InitialStateError {
    None,
    NonFiniteInput,
    InvalidObserverFrame,
    InvalidLocalDirection,
    SuperluminalLocalVelocity,
    NonFutureDirected,
    ConstraintViolation,
};

struct InitialStateResult {
    InitialStateError error = InitialStateError::None;
    std::optional<PhaseSpaceState> state;
    double measured_frequency =
        std::numeric_limits<double>::quiet_NaN();
    std::string message;

    explicit operator bool() const noexcept {
        return error == InitialStateError::None &&
               state.has_value();
    }
};

double observer_measured_frequency(
    const Covariant4& momentum,
    const Contravariant4& observer_velocity) noexcept;

InitialStateResult initialize_local_photon(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec3& local_direction,
    double affine = 0.0);

InitialStateResult initialize_local_timelike(
    const Metric& metric,
    const ObserverFrame& observer,
    const Vec3& local_velocity,
    double affine = 0.0);

} // namespace solar::relativity
