#include "kerr_separated_state.h"

#include "solar/relativity/hamiltonian.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

bool valid_direction(SeparatedDirection direction) noexcept {
    return direction == SeparatedDirection::Negative ||
           direction == SeparatedDirection::Locked ||
           direction == SeparatedDirection::Positive;
}

double direction_value(SeparatedDirection direction) {
    if (!valid_direction(direction)) {
        throw std::invalid_argument(
            "Kerr separated direction is not recognized");
    }
    return static_cast<double>(
        static_cast<int>(direction));
}

SeparatedDirection sign_direction(double value) {
    if (!std::isfinite(value) || value == 0.0) {
        throw std::domain_error(
            "Kerr separated direction has no finite sign");
    }
    return value < 0.0
               ? SeparatedDirection::Negative
               : SeparatedDirection::Positive;
}

bool finite_state(const PhaseSpaceState& state) noexcept {
    return std::isfinite(state.affine) &&
           state.x.v.all_finite() &&
           state.p.v.all_finite();
}

bool finite_mino_state(const KerrMinoState& state) noexcept {
    return std::all_of(
        state.begin(),
        state.end(),
        [](double value) {
            return std::isfinite(value);
        });
}

double normalized_radial_derivative(
    double derivative,
    double mass_M) noexcept {
    return std::fabs(derivative) /
           (mass_M * mass_M * mass_M);
}

double normalized_polar_derivative(
    double derivative,
    double mass_M) noexcept {
    return std::fabs(derivative) /
           (mass_M * mass_M);
}

SeparatedDirection classify_direction(
    double velocity,
    double potential,
    double potential_scale,
    double normalized_derivative,
    double physical_derivative,
    double potential_tolerance,
    double critical_tolerance,
    double integration_direction,
    bool symmetry_locked) {
    const double velocity_scale =
        std::sqrt(potential_scale);
    const double zero_velocity_threshold =
        std::sqrt(potential_tolerance) * velocity_scale;
    if (std::fabs(velocity) > zero_velocity_threshold) {
        return sign_direction(velocity);
    }

    if (std::fabs(potential) / potential_scale >
        potential_tolerance) {
        throw std::domain_error(
            "zero Kerr separated velocity is inconsistent "
            "with its potential");
    }
    if (symmetry_locked) {
        return SeparatedDirection::Locked;
    }
    if (normalized_derivative <= critical_tolerance) {
        throw KerrSeparatedCriticalInitialState(
            "Kerr separated initial state is near a "
            "critical double root");
    }
    return sign_direction(
        integration_direction * physical_derivative);
}

} // namespace

KerrSeparatedInitialState initialize_kerr_separated_state(
    const KerrBoyerLindquistMetric& metric,
    const PhaseSpaceState& initial,
    GeodesicKind kind,
    double normalized_potential_tolerance,
    double normalized_critical_tolerance,
    double integration_direction) {
    if (!is_valid_geodesic_kind(kind)) {
        throw std::invalid_argument(
            "Kerr separated geodesic kind is not recognized");
    }
    if (!std::isfinite(normalized_potential_tolerance) ||
        normalized_potential_tolerance <= 0.0 ||
        !std::isfinite(normalized_critical_tolerance) ||
        normalized_critical_tolerance <= 0.0 ||
        !std::isfinite(integration_direction) ||
        integration_direction == 0.0) {
        throw std::invalid_argument(
            "Kerr separated initialization tolerances and "
            "direction must be finite and non-zero");
    }
    if (!finite_state(initial)) {
        throw std::domain_error(
            "Kerr separated initial state is non-finite");
    }
    if (!metric.valid_point(initial.x)) {
        throw std::domain_error(
            "Kerr separated initial metric point is invalid");
    }

    const KerrConstants constants =
        evaluate_kerr_constants(metric, initial, kind);
    const double mu = std::cos(initial.x.v[2]);
    const KerrSeparatedPotentials potentials(
        metric.mass(), metric.spin_length(), constants);
    const KerrSeparatedPotentialValues values =
        potentials.evaluate(initial.x.v[1], mu);
    const PhaseSpaceDerivative hamiltonian_rhs =
        HamiltonGeodesicRhs(metric)(initial);
    const double radial_velocity =
        values.sigma * hamiltonian_rhs.dx.v[1];
    const double polar_velocity =
        -std::sin(initial.x.v[2]) *
        values.sigma *
        hamiltonian_rhs.dx.v[2];
    if (!std::isfinite(radial_velocity) ||
        !std::isfinite(polar_velocity)) {
        throw std::domain_error(
            "Kerr separated initial tangent is non-finite");
    }

    const double radial_residual =
        std::fabs(
            radial_velocity * radial_velocity -
            values.radial) /
        values.radial_scale;
    const double polar_residual =
        std::fabs(
            polar_velocity * polar_velocity -
            values.polar) /
        values.polar_scale;
    if (!std::isfinite(radial_residual) ||
        !std::isfinite(polar_residual) ||
        radial_residual > normalized_potential_tolerance ||
        polar_residual > normalized_potential_tolerance) {
        throw std::domain_error(
            "Kerr Hamiltonian state is inconsistent with "
            "the separated potentials");
    }

    const bool equatorial_lock =
        std::fabs(mu) <=
            std::sqrt(normalized_potential_tolerance) &&
        initial.p.v[2] == 0.0;
    const SeparatedDirection radial_direction =
        classify_direction(
            radial_velocity,
            values.radial,
            values.radial_scale,
            normalized_radial_derivative(
                values.radial_derivative, metric.mass()),
            values.radial_derivative,
            normalized_potential_tolerance,
            normalized_critical_tolerance,
            integration_direction,
            false);
    const SeparatedDirection polar_direction =
        classify_direction(
            polar_velocity,
            values.polar,
            values.polar_scale,
            normalized_polar_derivative(
                values.polar_derivative, metric.mass()),
            values.polar_derivative,
            normalized_potential_tolerance,
            normalized_critical_tolerance,
            integration_direction,
            equatorial_lock);

    return KerrSeparatedInitialState{
        KerrSeparatedState{
            KerrMinoState{{
                initial.x.v[0],
                initial.x.v[1],
                mu,
                initial.x.v[3],
                initial.affine,
            }},
            radial_direction,
            polar_direction,
        },
        constants,
    };
}

PhaseSpaceState reconstruct_kerr_phase_space(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& state) {
    if (!finite_mino_state(state.values) ||
        !valid_direction(state.radial_direction) ||
        !valid_direction(state.polar_direction)) {
        throw std::domain_error(
            "Kerr separated state cannot be reconstructed");
    }

    const double mu = state.values[kMu];
    if (std::fabs(mu) >= 1.0) {
        throw std::domain_error(
            "Kerr separated BL polar axis is unsupported");
    }
    const double theta = std::acos(mu);
    PhaseSpaceState result{
        state.values[kAffine],
        Contravariant4{
            Vec4{{
                state.values[kTime],
                state.values[kRadius],
                theta,
                state.values[kAzimuth],
            }}},
        Covariant4{
            Vec4{{-constants.E, 0.0, 0.0, constants.Lz}}},
    };
    if (!metric.valid_point(result.x)) {
        throw std::domain_error(
            "reconstructed Kerr BL point is invalid");
    }

    const KerrSeparatedPotentials potentials(
        metric.mass(), metric.spin_length(), constants);
    const KerrSeparatedPotentialValues values =
        potentials.evaluate(state.values[kRadius], mu);
    const double radial_sign =
        direction_value(state.radial_direction);
    const double polar_sign =
        direction_value(state.polar_direction);

    if (state.radial_direction == SeparatedDirection::Locked) {
        if (std::fabs(values.radial) >
            64.0 * std::numeric_limits<double>::epsilon() *
                values.radial_scale) {
            throw std::domain_error(
                "locked radial direction is not at a root");
        }
        result.p.v[1] = 0.0;
    } else {
        if (values.radial < 0.0 || values.delta == 0.0) {
            throw std::domain_error(
                "radial potential cannot be reconstructed");
        }
        result.p.v[1] =
            radial_sign * std::sqrt(values.radial) /
            values.delta;
    }

    if (state.polar_direction == SeparatedDirection::Locked) {
        if (std::fabs(values.polar) >
            64.0 * std::numeric_limits<double>::epsilon() *
                values.polar_scale) {
            throw std::domain_error(
                "locked polar direction is not at a root");
        }
        result.p.v[2] = 0.0;
    } else {
        const double sine = std::sqrt(1.0 - mu * mu);
        if (values.polar < 0.0 || sine == 0.0) {
            throw std::domain_error(
                "polar potential cannot be reconstructed");
        }
        result.p.v[2] =
            -polar_sign * std::sqrt(values.polar) / sine;
    }

    if (!finite_state(result)) {
        throw std::domain_error(
            "reconstructed Kerr phase-space state is non-finite");
    }
    return result;
}

} // namespace solar::relativity::detail
