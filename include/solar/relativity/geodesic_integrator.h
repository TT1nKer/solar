#pragma once

#include "solar/relativity/event_root.h"
#include "solar/relativity/hamiltonian.h"

#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

namespace solar::relativity {

using InvariantEvaluator =
    std::function<double(const PhaseSpaceState&)>;

struct GeodesicIntegrationConfig {
    GeodesicKind kind;
    numerics::Dopri5Config<8> dopri5;
    double initial_step;
    double min_step;
    double max_step;
    std::size_t max_rejections_per_step;
    std::size_t max_total_steps;
    double max_affine;
    double max_proper_time =
        std::numeric_limits<double>::infinity();
    double max_coordinate_time =
        std::numeric_limits<double>::infinity();
    double constraint_tolerance;
    bool monitor_energy = false;
    bool monitor_lz = false;
    InvariantEvaluator carter_evaluator;

    static GeodesicIntegrationConfig cpu_reference(
        GeodesicKind kind,
        double mass_scale,
        double initial_step,
        double max_step,
        double max_affine);
};

struct GeodesicIntegrationResult {
    PhaseSpaceState final_state;
    IntegrationDiagnostics diagnostics;
    std::optional<EventHit> event;
};

class GeodesicIntegrator {
public:
    explicit GeodesicIntegrator(const Metric& metric) noexcept;

    GeodesicIntegrationResult integrate(
        const PhaseSpaceState& initial,
        const GeodesicIntegrationConfig& config,
        const std::vector<GeodesicEvent>& events = {}) const;

private:
    const Metric* metric_;
};

} // namespace solar::relativity
