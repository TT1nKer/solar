#pragma once

#include "kerr_separated_potentials.h"
#include "kerr_separated_state.h"

#include "solar/relativity/kerr_separated.h"

#include <cstddef>
#include <optional>
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

using KerrTurningPhaseState = numerics::StateN<7>;

inline constexpr std::size_t kPhaseTime = 0;
inline constexpr std::size_t kPhaseRadius = 1;
inline constexpr std::size_t kPhaseRadialVelocity = 2;
inline constexpr std::size_t kPhaseMu = 3;
inline constexpr std::size_t kPhasePolarVelocity = 4;
inline constexpr std::size_t kPhaseAzimuth = 5;
inline constexpr std::size_t kPhaseAffine = 6;
// At a simple root, v^2 - V subtracts two roundoff-scale values.
// This remains bounded by the CPU geodesic constraint gate.
inline constexpr double kPhasePotentialDriftTolerance = 1.0e-10;

struct TurningPhaseCrossing {
    TurningStatus status;
    TurningCoordinate coordinate;
    double mino_parameter;
    double root_radius_M;
    std::string message;
};

struct TurningRelease {
    KerrSeparatedState state;
    KerrTurningPhaseState phase_state;
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
    const KerrTurningPhaseState& phase_state,
    double integration_direction,
    const KerrSeparatedConfig& config);

KerrTurningPhaseState initialize_kerr_turning_phase(
    const KerrSeparatedState& current,
    const KerrSeparatedPotentials& potentials);

KerrTurningPhaseState advance_kerr_turning_phase(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrTurningPhaseState& current,
    double mino_step,
    const KerrSeparatedConfig& config);

numerics::Dopri5StepResult<7>
attempt_kerr_turning_phase_step(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrTurningPhaseState& current,
    double current_mino,
    double attempted_step,
    const KerrSeparatedConfig& config);

KerrSeparatedState project_kerr_turning_phase(
    const KerrTurningPhaseState& phase,
    const KerrSeparatedState& direction_fallback);

std::optional<TurningPhaseCrossing>
locate_kerr_turning_phase_crossing(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentials& potentials,
    const numerics::Dopri5DenseOutput<7>& dense_output,
    double interval_end_mino,
    double normalized_potential_tolerance,
    double normalized_critical_derivative_tolerance);

} // namespace solar::relativity::detail
