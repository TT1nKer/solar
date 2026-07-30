#include "kerr_separated_turning.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace solar::relativity::detail {
namespace {

double velocity_for(
    TurningCoordinate coordinate,
    const KerrTurningPhaseState& state) noexcept {
    return state[
        coordinate == TurningCoordinate::Radial
            ? kPhaseRadialVelocity
            : kPhasePolarVelocity];
}

double potential_for(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentialValues& values) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? values.radial
               : values.polar;
}

double potential_scale_for(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentialValues& values) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? values.radial_scale
               : values.polar_scale;
}

double derivative_for(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentialValues& values) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? values.radial_derivative
               : values.polar_derivative;
}

double normalized_derivative(
    TurningCoordinate coordinate,
    double derivative,
    double mass_M) noexcept {
    const double scale =
        coordinate == TurningCoordinate::Radial
            ? mass_M * mass_M * mass_M
            : mass_M * mass_M;
    return std::fabs(derivative) / scale;
}

} // namespace

std::optional<TurningPhaseCrossing>
locate_kerr_turning_phase_crossing(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentials& potentials,
    const numerics::Dopri5DenseOutput<7>& dense_output,
    double interval_end_mino,
    double normalized_potential_tolerance,
    double normalized_critical_derivative_tolerance) {
    const double interval_start = dense_output.start();
    const double dense_end = dense_output.end();
    const double lower =
        std::min(interval_start, dense_end);
    const double upper =
        std::max(interval_start, dense_end);
    if (!std::isfinite(interval_end_mino) ||
        interval_end_mino < lower ||
        interval_end_mino > upper ||
        !std::isfinite(normalized_potential_tolerance) ||
        normalized_potential_tolerance <= 0.0 ||
        !std::isfinite(
            normalized_critical_derivative_tolerance) ||
        normalized_critical_derivative_tolerance <= 0.0) {
        return TurningPhaseCrossing{
            TurningStatus::Failed,
            coordinate,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            "turning phase crossing inputs are invalid",
        };
    }

    const auto start_state =
        dense_output.evaluate(interval_start);
    const auto end_state =
        dense_output.evaluate(interval_end_mino);
    double left_velocity =
        velocity_for(coordinate, start_state);
    double right_velocity =
        velocity_for(coordinate, end_state);
    if (!std::isfinite(left_velocity) ||
        !std::isfinite(right_velocity)) {
        return TurningPhaseCrossing{
            TurningStatus::Failed,
            coordinate,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            "turning phase velocity is non-finite",
        };
    }
    if (left_velocity != 0.0 &&
        right_velocity != 0.0 &&
        std::signbit(left_velocity) ==
            std::signbit(right_velocity)) {
        return std::nullopt;
    }

    double left = interval_start;
    double right = interval_end_mino;
    KerrTurningPhaseState root_state =
        std::fabs(left_velocity) <=
                std::fabs(right_velocity)
            ? start_state
            : end_state;
    double root_mino =
        std::fabs(left_velocity) <=
                std::fabs(right_velocity)
            ? left
            : right;
    for (std::size_t iteration = 0;
         iteration < 80 &&
         left_velocity != 0.0 &&
         right_velocity != 0.0;
         ++iteration) {
        if (std::nextafter(
                std::min(left, right),
                std::max(left, right)) ==
            std::max(left, right)) {
            break;
        }
        const double midpoint = 0.5 * (left + right);
        const auto midpoint_state =
            dense_output.evaluate(midpoint);
        const double midpoint_velocity =
            velocity_for(coordinate, midpoint_state);
        if (!std::isfinite(midpoint_velocity)) {
            return TurningPhaseCrossing{
                TurningStatus::Failed,
                coordinate,
                midpoint,
                std::numeric_limits<double>::quiet_NaN(),
                "turning phase midpoint velocity is non-finite",
            };
        }
        if (std::fabs(midpoint_velocity) <
            std::fabs(
                velocity_for(coordinate, root_state))) {
            root_state = midpoint_state;
            root_mino = midpoint;
        }
        if (midpoint_velocity == 0.0) {
            root_state = midpoint_state;
            root_mino = midpoint;
            break;
        }
        if (std::signbit(left_velocity) !=
            std::signbit(midpoint_velocity)) {
            right = midpoint;
            right_velocity = midpoint_velocity;
        } else {
            left = midpoint;
            left_velocity = midpoint_velocity;
        }
    }

    const auto root_values = potentials.evaluate(
        root_state[kPhaseRadius],
        root_state[kPhaseMu]);
    const double normalized_potential =
        std::fabs(
            potential_for(coordinate, root_values)) /
        potential_scale_for(coordinate, root_values);
    if (!std::isfinite(normalized_potential) ||
        normalized_potential >
            std::max(
                normalized_potential_tolerance,
                1.0e-10)) {
        return TurningPhaseCrossing{
            TurningStatus::Failed,
            coordinate,
            root_mino,
            root_state[kPhaseRadius],
            "turning phase crossing violates its potential",
        };
    }

    const double derivative =
        normalized_derivative(
            coordinate,
            derivative_for(coordinate, root_values),
            potentials.mass());
    const TurningStatus status =
        derivative <=
                normalized_critical_derivative_tolerance
            ? TurningStatus::NearCritical
            : TurningStatus::Simple;
    return TurningPhaseCrossing{
        status,
        coordinate,
        root_mino,
        root_state[kPhaseRadius],
        status == TurningStatus::Simple
            ? "simple turning phase crossing found"
            : "near-critical turning phase crossing found",
    };
}

} // namespace solar::relativity::detail
