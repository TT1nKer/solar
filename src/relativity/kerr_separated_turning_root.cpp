#include "kerr_separated_turning.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>

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

double potential_scale_for(
    TurningCoordinate coordinate,
    const KerrSeparatedPotentialValues& values) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? values.radial_scale
               : values.polar_scale;
}

KerrSeparatedPotentialValues evaluate_at(
    TurningCoordinate coordinate,
    double active_coordinate,
    double fixed_other_coordinate,
    const KerrSeparatedPotentials& potentials) {
    return coordinate == TurningCoordinate::Radial
               ? potentials.evaluate(
                     active_coordinate,
                     fixed_other_coordinate)
               : potentials.evaluate(
                     fixed_other_coordinate,
                     active_coordinate);
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

TurningRoot locate_kerr_turning_root(
    TurningCoordinate coordinate,
    double allowed_coordinate,
    double forbidden_coordinate,
    const KerrSeparatedPotentials& potentials,
    double fixed_other_coordinate,
    double normalized_root_tolerance,
    double normalized_potential_tolerance,
    double normalized_critical_derivative_tolerance) {
    if (!std::isfinite(allowed_coordinate) ||
        !std::isfinite(forbidden_coordinate) ||
        allowed_coordinate == forbidden_coordinate ||
        !std::isfinite(fixed_other_coordinate) ||
        !std::isfinite(normalized_root_tolerance) ||
        normalized_root_tolerance <= 0.0 ||
        !std::isfinite(normalized_potential_tolerance) ||
        normalized_potential_tolerance <= 0.0 ||
        !std::isfinite(
            normalized_critical_derivative_tolerance) ||
        normalized_critical_derivative_tolerance <= 0.0) {
        return TurningRoot{
            TurningStatus::Failed,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            0,
            "turning-root inputs are invalid",
        };
    }

    double allowed = allowed_coordinate;
    double forbidden = forbidden_coordinate;
    KerrSeparatedPotentialValues allowed_values;
    KerrSeparatedPotentialValues forbidden_values;
    try {
        allowed_values = evaluate_at(
            coordinate,
            allowed,
            fixed_other_coordinate,
            potentials);
        forbidden_values = evaluate_at(
            coordinate,
            forbidden,
            fixed_other_coordinate,
            potentials);
    } catch (const std::exception& error) {
        return TurningRoot{
            TurningStatus::Failed,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            0,
            std::string("turning-root evaluation failed: ") +
                error.what(),
        };
    }
    if (potential_for(coordinate, allowed_values) < 0.0 ||
        potential_for(coordinate, forbidden_values) >= 0.0) {
        return TurningRoot{
            TurningStatus::Failed,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            0,
            "turning root is not bracketed by allowed and forbidden points",
        };
    }

    const double coordinate_scale =
        coordinate == TurningCoordinate::Radial
            ? potentials.mass()
            : 1.0;
    const double machine_tolerance =
        64.0 * std::numeric_limits<double>::epsilon();
    const double internal_root_tolerance = std::min(
        normalized_root_tolerance, machine_tolerance);
    const double internal_potential_tolerance = std::min(
        normalized_potential_tolerance, machine_tolerance);
    std::size_t iterations = 0;
    for (; iterations < 100; ++iterations) {
        const double allowed_potential =
            potential_for(coordinate, allowed_values);
        const double allowed_potential_scale =
            potential_scale_for(coordinate, allowed_values);
        if (std::fabs(forbidden - allowed) /
                    coordinate_scale <=
                internal_root_tolerance &&
            std::fabs(allowed_potential) /
                    allowed_potential_scale <=
                internal_potential_tolerance) {
            break;
        }
        const double lower =
            std::min(allowed, forbidden);
        const double upper =
            std::max(allowed, forbidden);
        if (std::nextafter(lower, upper) == upper) {
            break;
        }
        const double forbidden_potential =
            potential_for(coordinate, forbidden_values);
        const double denominator =
            forbidden_potential - allowed_potential;
        double candidate =
            forbidden -
            forbidden_potential *
                (forbidden - allowed) / denominator;
        const double width = upper - lower;
        if (!std::isfinite(candidate) ||
            candidate <= lower + 0.1 * width ||
            candidate >= upper - 0.1 * width) {
            candidate = 0.5 * (allowed + forbidden);
        }

        KerrSeparatedPotentialValues candidate_values;
        try {
            candidate_values = evaluate_at(
                coordinate,
                candidate,
                fixed_other_coordinate,
                potentials);
        } catch (const std::exception& error) {
            return TurningRoot{
                TurningStatus::Failed,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                iterations,
                std::string(
                    "turning-root iteration failed: ") +
                    error.what(),
            };
        }
        if (potential_for(
                coordinate, candidate_values) >= 0.0) {
            allowed = candidate;
            allowed_values = candidate_values;
        } else {
            forbidden = candidate;
            forbidden_values = candidate_values;
        }
    }
    if (iterations == 100) {
        return TurningRoot{
            TurningStatus::Failed,
            allowed,
            std::numeric_limits<double>::quiet_NaN(),
            iterations,
            "turning-root iteration limit exceeded",
        };
    }

    double root_coordinate = allowed;
    KerrSeparatedPotentialValues root_values = allowed_values;
    double root_potential =
        potential_for(coordinate, root_values);
    const double forbidden_potential =
        potential_for(coordinate, forbidden_values);
    if (std::fabs(forbidden_potential) <
        std::fabs(root_potential)) {
        root_coordinate = forbidden;
        root_values = forbidden_values;
        root_potential = forbidden_potential;
    }
    const double denominator =
        forbidden_potential -
        potential_for(coordinate, allowed_values);
    const double secant =
        forbidden -
        forbidden_potential *
            (forbidden - allowed) / denominator;
    const double lower =
        std::min(allowed, forbidden);
    const double upper =
        std::max(allowed, forbidden);
    if (std::isfinite(secant) &&
        secant >= lower && secant <= upper) {
        try {
            const auto secant_values = evaluate_at(
                coordinate,
                secant,
                fixed_other_coordinate,
                potentials);
            const double secant_potential =
                potential_for(
                    coordinate, secant_values);
            if (std::fabs(secant_potential) <
                std::fabs(root_potential)) {
                root_coordinate = secant;
                root_values = secant_values;
                root_potential = secant_potential;
            }
        } catch (const std::exception&) {
            // The already-bracketed finite endpoints remain authoritative.
        }
    }
    const double derivative =
        derivative_for(coordinate, root_values);
    const double derivative_normalized =
        normalized_derivative(
            coordinate, derivative, potentials.mass());
    if (std::fabs(root_potential) /
                potential_scale_for(
                    coordinate, root_values) >
            normalized_potential_tolerance) {
        return TurningRoot{
            TurningStatus::Failed,
            root_coordinate,
            derivative_normalized,
            iterations,
            "located turning root does not meet the potential tolerance",
        };
    }
    const TurningStatus status =
        derivative_normalized <=
                normalized_critical_derivative_tolerance
            ? TurningStatus::NearCritical
            : TurningStatus::Simple;
    return TurningRoot{
        status,
        root_coordinate,
        derivative_normalized,
        iterations,
        status == TurningStatus::Simple
            ? "simple turning root found"
            : "near-critical turning root found",
    };
}

} // namespace solar::relativity::detail
