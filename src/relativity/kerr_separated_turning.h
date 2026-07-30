#pragma once

#include "kerr_separated_potentials.h"
#include "kerr_separated_state.h"

#include "solar/relativity/kerr_separated.h"

#include <cstddef>
#include <string>

namespace solar::relativity::detail {

enum class TurningCoordinate {
    Radial,
    Polar,
};

enum class TurningStatus {
    Simple,
    NearCritical,
    Failed,
};

struct TurningRoot {
    TurningStatus status;
    double coordinate;
    double normalized_derivative;
    std::size_t iterations;
    std::string message;
};

struct TurningRelease {
    KerrSeparatedState state;
    double mino_step;
    double root_radius_M;
};

TurningRoot locate_kerr_turning_root(
    TurningCoordinate coordinate,
    double allowed_coordinate,
    double forbidden_coordinate,
    const KerrSeparatedPotentials& potentials,
    double fixed_other_coordinate,
    double normalized_root_tolerance,
    double normalized_potential_tolerance,
    double normalized_critical_derivative_tolerance);

TurningRelease release_kerr_turning_point(
    TurningCoordinate coordinate,
    const TurningRoot& root,
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& current,
    double integration_direction,
    const KerrSeparatedConfig& config);

} // namespace solar::relativity::detail
