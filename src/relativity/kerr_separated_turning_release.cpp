#include "kerr_separated_turning.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace solar::relativity::detail {
namespace {

double potential_for(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentialValues& values) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? values.radial
               : values.polar;
}

double derivative_for(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentialValues& values) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? values.radial_derivative
               : values.polar_derivative;
}

double direction_sign(SeparatedDirection direction) {
    if (direction == SeparatedDirection::Negative) {
        return -1.0;
    }
    if (direction == SeparatedDirection::Positive) {
        return 1.0;
    }
    if (direction == SeparatedDirection::Locked) {
        return 0.0;
    }
    throw std::invalid_argument(
        "turning direction is not recognized");
}

SeparatedDirection direction_from_velocity(double velocity) {
    if (!std::isfinite(velocity) || velocity == 0.0) {
        throw std::domain_error(
            "turning velocity has no finite direction");
    }
    return velocity < 0.0
               ? SeparatedDirection::Negative
               : SeparatedDirection::Positive;
}

double active_phase_coordinate(
    TurningCoordinate coordinate,
    const KerrTurningPhaseState& state) noexcept {
    return state[
        coordinate == TurningCoordinate::Radial
            ? kPhaseRadius
            : kPhaseMu];
}

double active_phase_velocity(
    TurningCoordinate coordinate,
    const KerrTurningPhaseState& state) noexcept {
    return state[
        coordinate == TurningCoordinate::Radial
            ? kPhaseRadialVelocity
            : kPhasePolarVelocity];
}

void require_phase_constraint(
    const KerrTurningPhaseState& phase,
    const KerrSeparatedPotentials& potentials,
    double normalized_tolerance) {
    const auto values = potentials.evaluate(
        phase[kPhaseRadius], phase[kPhaseMu]);
    const double radial_residual =
        std::fabs(
            phase[kPhaseRadialVelocity] *
                phase[kPhaseRadialVelocity] -
            values.radial) /
        values.radial_scale;
    const double polar_residual =
        std::fabs(
            phase[kPhasePolarVelocity] *
                phase[kPhasePolarVelocity] -
            values.polar) /
        values.polar_scale;
    if (!std::isfinite(radial_residual) ||
        !std::isfinite(polar_residual) ||
        radial_residual > normalized_tolerance ||
        polar_residual > normalized_tolerance) {
        throw std::domain_error(
            "turning phase drift exceeds the potential tolerance");
    }
}

} // namespace

TurningRelease release_kerr_turning_point(
    TurningCoordinate coordinate,
    const TurningRoot& root,
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& current,
    const KerrTurningPhaseState& phase_state,
    double integration_direction,
    const KerrSeparatedConfig& config) {
    if (root.status != TurningStatus::Simple ||
        !std::isfinite(root.coordinate) ||
        !std::isfinite(integration_direction) ||
        integration_direction == 0.0) {
        throw std::invalid_argument(
            "simple turning root and direction are required");
    }

    const KerrSeparatedPotentials potentials(
        metric.mass(), metric.spin_length(), constants);
    require_phase_constraint(
        phase_state,
        potentials,
        std::max(
            config.potential_tolerance,
            kPhasePotentialDriftTolerance));
    const double fixed_other_coordinate =
        phase_state[
            coordinate == TurningCoordinate::Radial
                ? kPhaseMu
                : kPhaseRadius];
    const auto root_values =
        coordinate == TurningCoordinate::Radial
            ? potentials.evaluate(
                  root.coordinate,
                  fixed_other_coordinate)
            : potentials.evaluate(
                  fixed_other_coordinate,
                  root.coordinate);
    const double physical_derivative =
        derivative_for(coordinate, root_values);
    if (!std::isfinite(physical_derivative) ||
        physical_derivative == 0.0) {
        throw std::domain_error(
            "simple turning root has no release derivative");
    }

    const SeparatedDirection incoming_direction =
        coordinate == TurningCoordinate::Radial
            ? current.radial_direction
            : current.polar_direction;
    const double incoming_sign =
        direction_sign(incoming_direction);
    const double phase_velocity =
        active_phase_velocity(coordinate, phase_state);
    double approach_step = 0.0;
    if (phase_velocity * incoming_sign > 0.0) {
        approach_step =
            -2.0 * phase_velocity / physical_derivative;
        if (!std::isfinite(approach_step) ||
            approach_step * integration_direction < 0.0) {
            throw std::domain_error(
                "turning phase approach has the wrong direction");
        }
    }

    const double coordinate_tolerance =
        config.root_tolerance *
        (coordinate == TurningCoordinate::Radial
             ? metric.mass()
             : 1.0);
    const double release_magnitude = std::min(
        config.max_mino_step,
        std::max(
            config.min_mino_step,
            2.0 *
                std::sqrt(
                    coordinate_tolerance /
                    std::fabs(physical_derivative))));
    const double transition_step =
        approach_step +
        std::copysign(
            release_magnitude, integration_direction);
    const KerrTurningPhaseState transitioned =
        advance_kerr_turning_phase(
            metric,
            constants,
            potentials,
            phase_state,
            transition_step,
            config);
    require_phase_constraint(
        transitioned,
        potentials,
        std::max(
            config.potential_tolerance,
            kPhasePotentialDriftTolerance));

    KerrSeparatedState released =
        project_kerr_turning_phase(transitioned, current);
    const auto released_values = potentials.evaluate(
        released.values[kRadius], released.values[kMu]);
    if (potential_for(coordinate, released_values) < 0.0) {
        throw std::domain_error(
            "turning release did not enter the allowed region");
    }
    const SeparatedDirection expected_direction =
        direction_from_velocity(
            physical_derivative * integration_direction);
    const SeparatedDirection actual_direction =
        coordinate == TurningCoordinate::Radial
            ? released.radial_direction
            : released.polar_direction;
    if (actual_direction != expected_direction) {
        throw std::domain_error(
            "turning transition did not cross the root");
    }

    const double root_radius =
        coordinate == TurningCoordinate::Radial
            ? root.coordinate
            : std::min(
                  active_phase_coordinate(
                      TurningCoordinate::Radial,
                      phase_state),
                  released.values[kRadius]);
    return TurningRelease{
        std::move(released),
        transitioned,
        transition_step,
        root_radius,
    };
}

} // namespace solar::relativity::detail
