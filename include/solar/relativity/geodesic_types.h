#pragma once

#include "solar/relativity/types.h"

#include <cstddef>
#include <functional>
#include <limits>
#include <string>

namespace solar::relativity {

enum class EventDirection {
    Any,
    Increasing,
    Decreasing,
};

constexpr bool is_valid_event_direction(
    EventDirection direction) noexcept {
    return direction == EventDirection::Any ||
           direction == EventDirection::Increasing ||
           direction == EventDirection::Decreasing;
}

enum class TerminationReason {
    HorizonCrossing,
    InteriorCutoff,
    Escaped,
    DiskSurfaceHit,
    MaterialSurfaceHit,
    RadialTurningPoint,
    PolarTurningPoint,
    MaxAffine,
    MaxProperTime,
    MaxCoordinateTime,
    MaxSteps,
    StepUnderflow,
    InvalidMetricPoint,
    NonFiniteState,
    ConstraintViolation,
    EventRootFailure,
    UserEvent,
    NearCriticalOrbit,
};

struct IntegrationDiagnostics {
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    double min_step = std::numeric_limits<double>::quiet_NaN();
    double max_step = std::numeric_limits<double>::quiet_NaN();
    double max_constraint_error = 0.0;
    double max_energy_rel_error =
        std::numeric_limits<double>::quiet_NaN();
    double max_lz_rel_error =
        std::numeric_limits<double>::quiet_NaN();
    double max_carter_rel_error =
        std::numeric_limits<double>::quiet_NaN();
    TerminationReason reason = TerminationReason::NonFiniteState;
    std::string message;
    double max_carter_abs_error =
        std::numeric_limits<double>::quiet_NaN();
};

using EventFunction =
    std::function<double(const PhaseSpaceState&)>;

struct GeodesicEvent {
    std::string name;
    EventFunction function;
    EventDirection direction = EventDirection::Any;
    TerminationReason reason = TerminationReason::UserEvent;
    double root_tolerance = 1.0e-10;
};

struct EventHit {
    std::size_t event_index;
    double affine;
    PhaseSpaceState state;
    double value;
    std::size_t root_iterations;
};

} // namespace solar::relativity
