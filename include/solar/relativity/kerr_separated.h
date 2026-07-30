#pragma once

#include "solar/numerics/dopri5.h"
#include "solar/relativity/geodesic_types.h"
#include "solar/relativity/kerr_constants.h"

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace solar::relativity {

struct KerrSeparatedConfig {
    GeodesicKind kind;
    numerics::Dopri5Config<5> dopri5;
    double initial_mino_step;
    double min_mino_step;
    double max_mino_step;
    std::size_t max_rejections_per_step;
    std::size_t max_total_steps;
    double max_affine;
    double max_coordinate_time;
    double potential_tolerance;
    double root_tolerance;
    double critical_derivative_tolerance;
    double polar_axis_tolerance;

    static KerrSeparatedConfig cpu_reference(
        GeodesicKind kind,
        double mass_scale,
        double initial_mino_step,
        double max_mino_step,
        double max_affine);
};

struct KerrSeparatedDiagnostics {
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t radial_turns = 0;
    std::size_t polar_turns = 0;
    double min_mino_step =
        std::numeric_limits<double>::quiet_NaN();
    double max_mino_step =
        std::numeric_limits<double>::quiet_NaN();
    double min_radius_M =
        std::numeric_limits<double>::quiet_NaN();
    double azimuthal_advance = 0.0;
    double winding = 0.0;
    double max_radial_residual = 0.0;
    double max_polar_residual = 0.0;
    double max_constraint_error = 0.0;
    double max_carter_rel_error = 0.0;
    TerminationReason reason = TerminationReason::NonFiniteState;
    std::string message;
};

struct KerrSeparatedIntegrationResult {
    PhaseSpaceState final_state;
    KerrConstants constants{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
    };
    KerrSeparatedDiagnostics diagnostics;
    std::optional<EventHit> event;
};

class KerrSeparatedIntegrator {
public:
    explicit KerrSeparatedIntegrator(
        const KerrBoyerLindquistMetric& metric) noexcept;

    KerrSeparatedIntegrationResult integrate(
        const PhaseSpaceState& initial,
        const KerrSeparatedConfig& config,
        const std::vector<GeodesicEvent>& events = {}) const;

private:
    const KerrBoyerLindquistMetric* metric_;
};

} // namespace solar::relativity
