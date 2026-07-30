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

SeparatedDirection direction_from(double value) {
    if (!std::isfinite(value) || value == 0.0) {
        throw std::domain_error(
            "turning release direction has no sign");
    }
    return value < 0.0
               ? SeparatedDirection::Negative
               : SeparatedDirection::Positive;
}

double signed_sqrt(
    SeparatedDirection direction,
    double potential) {
    if (direction == SeparatedDirection::Locked) {
        return 0.0;
    }
    if (potential < 0.0 || !std::isfinite(potential)) {
        throw std::domain_error(
            "turning release encountered forbidden potential");
    }
    return static_cast<double>(
               static_cast<int>(direction)) *
           std::sqrt(potential);
}

void set_active_coordinate(
    TurningCoordinate coordinate,
    KerrMinoState& state,
    double value) noexcept {
    state[coordinate == TurningCoordinate::Radial
              ? kRadius
              : kMu] = value;
}

double active_coordinate(
    TurningCoordinate coordinate,
    const KerrMinoState& state) noexcept {
    return state[coordinate == TurningCoordinate::Radial
                     ? kRadius
                     : kMu];
}

double other_coordinate(
    TurningCoordinate coordinate,
    const KerrMinoState& state) noexcept {
    return state[coordinate == TurningCoordinate::Radial
                     ? kMu
                     : kRadius];
}

void set_other_coordinate(
    TurningCoordinate coordinate,
    KerrMinoState& state,
    double value) noexcept {
    state[coordinate == TurningCoordinate::Radial
              ? kMu
              : kRadius] = value;
}

SeparatedDirection other_direction(
    TurningCoordinate coordinate,
    const KerrSeparatedState& state) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? state.polar_direction
               : state.radial_direction;
}

SeparatedDirection active_direction(
    TurningCoordinate coordinate,
    const KerrSeparatedState& state) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? state.radial_direction
               : state.polar_direction;
}

void set_outgoing_direction(
    TurningCoordinate coordinate,
    KerrSeparatedState& state,
    SeparatedDirection direction) noexcept {
    if (coordinate == TurningCoordinate::Radial) {
        state.radial_direction = direction;
    } else {
        state.polar_direction = direction;
    }
}

} // namespace

TurningRelease release_kerr_turning_point(
    TurningCoordinate coordinate,
    const TurningRoot& root,
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& current,
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
    const auto current_values = potentials.evaluate(
        current.values[kRadius], current.values[kMu]);
    KerrMinoState root_state = current.values;
    set_active_coordinate(
        coordinate, root_state, root.coordinate);
    const auto root_values = potentials.evaluate(
        root_state[kRadius], root_state[kMu]);
    const double physical_derivative =
        derivative_for(coordinate, root_values);
    if (!std::isfinite(physical_derivative) ||
        physical_derivative == 0.0) {
        throw std::domain_error(
            "simple turning root has no release derivative");
    }
    const double current_active_potential =
        potential_for(coordinate, current_values);
    const double incoming_velocity =
        signed_sqrt(
            active_direction(coordinate, current),
            current_active_potential);
    double approach_step =
        -2.0 * incoming_velocity / physical_derivative;
    if (!std::isfinite(approach_step)) {
        throw std::domain_error(
            "turning approach step is non-finite");
    }
    if (approach_step * integration_direction < 0.0) {
        if (std::fabs(approach_step) <=
            config.min_mino_step) {
            approach_step = 0.0;
        } else {
            throw std::domain_error(
                "turning approach has the wrong direction");
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
    const double release_step =
        std::copysign(
            release_magnitude, integration_direction);
    const double transition_step =
        approach_step + release_step;
    const SeparatedDirection outgoing =
        direction_from(
            physical_derivative * release_step);

    KerrMinoState midpoint = current.values;
    const double half_step = 0.5 * transition_step;
    set_active_coordinate(
        coordinate,
        midpoint,
        active_coordinate(
            coordinate, current.values) +
            incoming_velocity * half_step +
            0.25 * physical_derivative *
                half_step * half_step);

    const double root_other =
        other_coordinate(coordinate, current.values);
    const double other_current_potential =
        coordinate == TurningCoordinate::Radial
            ? current_values.polar
            : current_values.radial;
    const double other_current_velocity =
        signed_sqrt(
            other_direction(coordinate, current),
            other_current_potential);
    set_other_coordinate(
        coordinate,
        midpoint,
        root_other + half_step * other_current_velocity);
    const auto midpoint_values = potentials.evaluate(
        midpoint[kRadius], midpoint[kMu]);
    const double other_midpoint_potential =
        coordinate == TurningCoordinate::Radial
            ? midpoint_values.polar
            : midpoint_values.radial;
    const double other_midpoint_velocity =
        signed_sqrt(
            other_direction(coordinate, current),
            other_midpoint_potential);

    KerrSeparatedState released = current;
    released.values = current.values;
    set_active_coordinate(
        coordinate,
        released.values,
        root.coordinate +
            0.25 * physical_derivative *
                release_step * release_step);
    set_other_coordinate(
        coordinate,
        released.values,
        root_other +
            transition_step * other_midpoint_velocity);
    set_outgoing_direction(
        coordinate, released, outgoing);

    const double radius = midpoint[kRadius];
    const double mu = midpoint[kMu];
    const double radius_sq = radius * radius;
    const double spin = metric.spin_length();
    const double spin_sq = spin * spin;
    const double radial_momentum =
        constants.E * (radius_sq + spin_sq) -
        spin * constants.Lz;
    const double one_minus_mu_sq = 1.0 - mu * mu;
    const double time_derivative =
        (radius_sq + spin_sq) * radial_momentum /
            midpoint_values.delta +
        spin *
            (constants.Lz -
             spin * constants.E * one_minus_mu_sq);
    const double azimuth_derivative =
        spin * radial_momentum /
            midpoint_values.delta +
        constants.Lz / one_minus_mu_sq -
        spin * constants.E;
    released.values[kTime] =
        current.values[kTime] +
        transition_step * time_derivative;
    released.values[kAzimuth] =
        current.values[kAzimuth] +
        transition_step * azimuth_derivative;
    released.values[kAffine] =
        current.values[kAffine] +
        transition_step * midpoint_values.sigma;

    const auto released_values = potentials.evaluate(
        released.values[kRadius], released.values[kMu]);
    if (potential_for(coordinate, released_values) < 0.0) {
        throw std::domain_error(
            "turning release did not enter the allowed region");
    }
    const double root_radius =
        coordinate == TurningCoordinate::Radial
            ? root.coordinate
            : std::min(
                  current.values[kRadius],
                  released.values[kRadius]);
    return TurningRelease{
        std::move(released),
        transition_step,
        root_radius,
    };
}

} // namespace solar::relativity::detail
