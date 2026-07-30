#include "kerr_separated_turning.h"

#include "solar/numerics/dopri5.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

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

} // namespace solar::relativity::detail
