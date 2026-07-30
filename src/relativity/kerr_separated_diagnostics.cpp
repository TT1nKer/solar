#include "kerr_separated_diagnostics.h"

#include "kerr_separated_potentials.h"

#include "solar/relativity/hamiltonian.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace solar::relativity::detail {
namespace {

constexpr double kTwoPi =
    6.283185307179586476925286766559005768;

void require_finite(double value, const char* message) {
    if (!std::isfinite(value)) {
        throw std::domain_error(message);
    }
}

} // namespace

KerrSeparatedDiagnosticTracker::
    KerrSeparatedDiagnosticTracker(
        const KerrBoyerLindquistMetric& metric,
        GeodesicKind kind,
        KerrConstants constants,
        const PhaseSpaceState& initial) noexcept
    : metric_(&metric),
      kind_(kind),
      constants_(constants),
      initial_azimuth_(initial.x.v[3]),
      carter_denominator_(
          std::max(
              std::fabs(constants.Q),
              metric.mass() * metric.mass() *
                  1.0e-14)) {}

void KerrSeparatedDiagnosticTracker::initialize(
    const KerrSeparatedState& separated,
    const PhaseSpaceState& public_state,
    KerrSeparatedDiagnostics& diagnostics) const {
    diagnostics.min_radius_M = public_state.x.v[1];
    diagnostics.azimuthal_advance =
        public_state.x.v[3] - initial_azimuth_;
    diagnostics.winding =
        diagnostics.azimuthal_advance / kTwoPi;
    update_invariants(
        separated, public_state, diagnostics);
}

void KerrSeparatedDiagnosticTracker::accept(
    const KerrSeparatedState& separated,
    const PhaseSpaceState& public_state,
    double mino_step_magnitude,
    KerrSeparatedDiagnostics& diagnostics,
    double additional_minimum_radius_M) const {
    if (!std::isfinite(mino_step_magnitude) ||
        mino_step_magnitude <= 0.0) {
        throw std::domain_error(
            "accepted Mino step magnitude is invalid");
    }

    KerrSeparatedDiagnostics updated = diagnostics;
    update_invariants(separated, public_state, updated);
    ++updated.accepted_steps;
    updated.min_mino_step =
        std::isnan(updated.min_mino_step)
            ? mino_step_magnitude
            : std::min(
                  updated.min_mino_step,
                  mino_step_magnitude);
    updated.max_mino_step =
        std::isnan(updated.max_mino_step)
            ? mino_step_magnitude
            : std::max(
                  updated.max_mino_step,
                  mino_step_magnitude);
    updated.min_radius_M = std::min(
        updated.min_radius_M,
        public_state.x.v[1]);
    if (std::isfinite(additional_minimum_radius_M)) {
        updated.min_radius_M = std::min(
            updated.min_radius_M,
            additional_minimum_radius_M);
    }
    updated.azimuthal_advance =
        public_state.x.v[3] - initial_azimuth_;
    updated.winding =
        updated.azimuthal_advance / kTwoPi;
    diagnostics = std::move(updated);
}

void KerrSeparatedDiagnosticTracker::update_invariants(
    const KerrSeparatedState& separated,
    const PhaseSpaceState& public_state,
    KerrSeparatedDiagnostics& diagnostics) const {
    const double constraint =
        hamiltonian_constraint_error(
            *metric_, public_state, kind_);
    const KerrConstants measured =
        evaluate_kerr_constants(
            *metric_, public_state, kind_);
    const double carter_relative_error =
        std::fabs(measured.Q - constants_.Q) /
        carter_denominator_;

    const KerrSeparatedPotentials potentials(
        metric_->mass(),
        metric_->spin_length(),
        constants_);
    const auto values = potentials.evaluate(
        separated.values[kRadius],
        separated.values[kMu]);
    const PhaseSpaceDerivative tangent =
        HamiltonGeodesicRhs(*metric_)(public_state);
    const double radial_velocity =
        values.sigma * tangent.dx.v[1];
    const double polar_velocity =
        -std::sin(public_state.x.v[2]) *
        values.sigma * tangent.dx.v[2];
    const double radial_velocity_sq =
        radial_velocity * radial_velocity;
    const double polar_velocity_sq =
        polar_velocity * polar_velocity;
    const double radial_scale = std::max(
        values.radial_scale,
        std::fabs(radial_velocity_sq));
    const double polar_scale = std::max(
        values.polar_scale,
        std::fabs(polar_velocity_sq));
    const double radial_residual =
        std::fabs(radial_velocity_sq - values.radial) /
        radial_scale;
    const double polar_residual =
        std::fabs(polar_velocity_sq - values.polar) /
        polar_scale;

    require_finite(
        constraint,
        "Hamiltonian diagnostic is non-finite");
    require_finite(
        carter_relative_error,
        "Carter diagnostic is non-finite");
    require_finite(
        radial_residual,
        "radial potential diagnostic is non-finite");
    require_finite(
        polar_residual,
        "polar potential diagnostic is non-finite");

    diagnostics.max_constraint_error = std::max(
        diagnostics.max_constraint_error, constraint);
    diagnostics.max_carter_rel_error = std::max(
        diagnostics.max_carter_rel_error,
        carter_relative_error);
    diagnostics.max_radial_residual = std::max(
        diagnostics.max_radial_residual,
        radial_residual);
    diagnostics.max_polar_residual = std::max(
        diagnostics.max_polar_residual,
        polar_residual);
}

} // namespace solar::relativity::detail
