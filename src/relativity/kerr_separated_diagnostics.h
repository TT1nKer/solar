#pragma once

#include "kerr_separated_state.h"

#include "solar/relativity/kerr_separated.h"

namespace solar::relativity::detail {

class KerrSeparatedDiagnosticTracker {
public:
    KerrSeparatedDiagnosticTracker(
        const KerrBoyerLindquistMetric& metric,
        GeodesicKind kind,
        KerrConstants constants,
        const PhaseSpaceState& initial) noexcept;

    void initialize(
        const KerrSeparatedState& separated,
        const PhaseSpaceState& public_state,
        KerrSeparatedDiagnostics& diagnostics) const;

    void accept(
        const KerrSeparatedState& separated,
        const PhaseSpaceState& public_state,
        double mino_step_magnitude,
        KerrSeparatedDiagnostics& diagnostics,
        double additional_minimum_radius_M =
            std::numeric_limits<double>::quiet_NaN()) const;

private:
    void update_invariants(
        const KerrSeparatedState& separated,
        const PhaseSpaceState& public_state,
        KerrSeparatedDiagnostics& diagnostics) const;

    const KerrBoyerLindquistMetric* metric_;
    GeodesicKind kind_;
    KerrConstants constants_;
    double initial_azimuth_;
    double carter_denominator_;
};

} // namespace solar::relativity::detail
