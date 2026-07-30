#include "kerr_separated_turning.h"

#include "solar/numerics/dopri5.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

double initial_velocity(
    SeparatedDirection direction,
    double potential) {
    if (direction == SeparatedDirection::Locked) {
        return 0.0;
    }
    if (!std::isfinite(potential) || potential < 0.0) {
        throw std::domain_error(
            "turning phase starts outside the allowed potential");
    }
    return direction_sign(direction) *
           std::sqrt(potential);
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

numerics::Dopri5Config<7> phase_dopri_config(
    const KerrBoyerLindquistMetric& metric,
    const KerrSeparatedConfig& config) {
    const double mass = metric.mass();
    numerics::Dopri5Config<7> result{};
    result.absolute_tolerance = {{
        1.0e-15 * mass,
        1.0e-15 * mass,
        1.0e-14 * mass * mass,
        1.0e-16,
        1.0e-15 * mass,
        1.0e-16,
        1.0e-15 * mass,
    }};
    result.relative_tolerance =
        std::min(
            config.dopri5.relative_tolerance,
            2.0e-15);
    result.safety = config.dopri5.safety;
    result.min_factor = config.dopri5.min_factor;
    result.max_factor = config.dopri5.max_factor;
    result.error_norm = config.dopri5.error_norm;
    return result;
}

KerrTurningPhaseState phase_derivative(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrTurningPhaseState& state) {
    const double radius = state[kPhaseRadius];
    const double mu = state[kPhaseMu];
    if (!std::isfinite(radius) ||
        !std::isfinite(mu) ||
        std::fabs(mu) >= 1.0) {
        throw std::domain_error(
            "turning phase left the BL domain");
    }

    const auto values = potentials.evaluate(radius, mu);
    const double radius_sq = radius * radius;
    const double spin = metric.spin_length();
    const double spin_sq = spin * spin;
    const double one_minus_mu_sq = 1.0 - mu * mu;
    const double radial_momentum =
        constants.E * (radius_sq + spin_sq) -
        spin * constants.Lz;

    KerrTurningPhaseState derivative{};
    derivative[kPhaseTime] =
        (radius_sq + spin_sq) * radial_momentum /
            values.delta +
        spin *
            (constants.Lz -
             spin * constants.E * one_minus_mu_sq);
    derivative[kPhaseRadius] =
        state[kPhaseRadialVelocity];
    derivative[kPhaseRadialVelocity] =
        0.5 * values.radial_derivative;
    derivative[kPhaseMu] =
        state[kPhasePolarVelocity];
    derivative[kPhasePolarVelocity] =
        0.5 * values.polar_derivative;
    derivative[kPhaseAzimuth] =
        spin * radial_momentum / values.delta +
        constants.Lz / one_minus_mu_sq -
        spin * constants.E;
    derivative[kPhaseAffine] = values.sigma;
    return derivative;
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

KerrTurningPhaseState initialize_kerr_turning_phase(
    const KerrSeparatedState& current,
    const KerrSeparatedPotentials& potentials) {
    const auto values = potentials.evaluate(
        current.values[kRadius], current.values[kMu]);
    return KerrTurningPhaseState{{
        current.values[kTime],
        current.values[kRadius],
        initial_velocity(
            current.radial_direction, values.radial),
        current.values[kMu],
        initial_velocity(
            current.polar_direction, values.polar),
        current.values[kAzimuth],
        current.values[kAffine],
    }};
}

KerrTurningPhaseState advance_kerr_turning_phase(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrTurningPhaseState& current,
    double mino_step,
    const KerrSeparatedConfig& config) {
    if (!std::isfinite(mino_step) || mino_step == 0.0) {
        throw std::invalid_argument(
            "turning phase step must be finite and non-zero");
    }

    KerrTurningPhaseState state = current;
    double elapsed = 0.0;
    double proposed_step = mino_step;
    for (std::size_t attempt = 0;
         attempt < 512;
         ++attempt) {
        const double remaining = mino_step - elapsed;
        if (remaining == 0.0) {
            return state;
        }
        const double step = std::copysign(
            std::min(
                std::fabs(proposed_step),
                std::fabs(remaining)),
            mino_step);
        const auto result =
            attempt_kerr_turning_phase_step(
                metric,
                constants,
                potentials,
                state,
                elapsed,
                step,
                config);
        if (result.status !=
            numerics::Dopri5StepResult<7>::Status::Completed) {
            throw std::domain_error(
                "turning phase derivative is non-finite");
        }
        proposed_step = result.next_step;
        if (!result.accepted) {
            continue;
        }
        state = result.state;
        elapsed += step;
        if (elapsed == mino_step ||
            std::fabs(mino_step - elapsed) <=
                8.0 * std::numeric_limits<double>::epsilon() *
                    std::fabs(mino_step)) {
            return state;
        }
    }
    throw std::domain_error(
        "turning phase exceeded its adaptive-step budget");
}

numerics::Dopri5StepResult<7>
attempt_kerr_turning_phase_step(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrTurningPhaseState& current,
    double current_mino,
    double attempted_step,
    const KerrSeparatedConfig& config) {
    const auto rhs =
        [&](double, const KerrTurningPhaseState& state) {
            return phase_derivative(
                metric, constants, potentials, state);
        };
    return numerics::dopri5_step(
        current,
        current_mino,
        attempted_step,
        rhs,
        phase_dopri_config(metric, config));
}

KerrSeparatedState project_kerr_turning_phase(
    const KerrTurningPhaseState& phase,
    const KerrSeparatedState& direction_fallback) {
    return KerrSeparatedState{
        KerrMinoState{{
            phase[kPhaseTime],
            phase[kPhaseRadius],
            phase[kPhaseMu],
            phase[kPhaseAzimuth],
            phase[kPhaseAffine],
        }},
        direction_fallback.radial_direction ==
                SeparatedDirection::Locked
            ? SeparatedDirection::Locked
            : direction_from_velocity(
                  phase[kPhaseRadialVelocity]),
        direction_fallback.polar_direction ==
                SeparatedDirection::Locked
            ? SeparatedDirection::Locked
            : direction_from_velocity(
                  phase[kPhasePolarVelocity]),
    };
}

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
            1.0e-10));
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
            1.0e-10));

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
