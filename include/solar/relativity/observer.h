#pragma once

#include "solar/relativity/metric.h"

#include <array>
#include <optional>
#include <string>

namespace solar::relativity {

struct Tetrad {
    // basis[a] stores e_(a)^mu.
    std::array<Contravariant4, 4> basis;
};

struct ObserverFrame {
    Contravariant4 x;
    Tetrad tetrad;
};

enum class ObserverError {
    None,
    NonFiniteInput,
    InvalidMetricPoint,
    FourVelocityNotUnitTimelike,
    DegenerateSpatialSeed,
    StaticWorldlineNotTimelike,
    CircularWorldlineNotTimelike,
    TetradValidationFailure,
};

struct ObserverResult {
    ObserverError error = ObserverError::None;
    std::optional<ObserverFrame> frame;
    std::string message;

    explicit operator bool() const noexcept {
        return error == ObserverError::None &&
               frame.has_value();
    }
};

struct LookAtAttitude {
    Contravariant4 look_direction;
    Contravariant4 up_reference;
};

double tetrad_orthonormality_error(
    const Metric& metric,
    const ObserverFrame& observer);

Contravariant4 tetrad_to_coordinate(
    const Tetrad& tetrad,
    const Vec4& local_components) noexcept;

Vec4 coordinate_to_tetrad(
    const Mat4& covariant,
    const Tetrad& tetrad,
    const Contravariant4& coordinate_vector) noexcept;

ObserverResult make_static_observer(
    const Metric& metric,
    const Contravariant4& x);

ObserverResult make_arbitrary_observer(
    const Metric& metric,
    const Contravariant4& x,
    const Contravariant4& four_velocity,
    const std::array<Contravariant4, 3>& spatial_seeds);

ObserverResult make_look_at_observer(
    const Metric& metric,
    const Contravariant4& x,
    const Contravariant4& four_velocity,
    const LookAtAttitude& attitude);

} // namespace solar::relativity
