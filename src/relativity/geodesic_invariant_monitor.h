#pragma once

#include "solar/relativity/geodesic_integrator.h"

#include <string>

namespace solar::relativity::detail {

class GeodesicInvariantMonitor {
public:
    GeodesicInvariantMonitor(
        const PhaseSpaceState& initial,
        const GeodesicIntegrationConfig& config);

    bool initialize(
        IntegrationDiagnostics& diagnostics,
        std::string& failure_message);

    bool update(
        const PhaseSpaceState& state,
        IntegrationDiagnostics& diagnostics,
        std::string& failure_message) const;

private:
    bool evaluate_carter(
        const PhaseSpaceState& state,
        double& value,
        std::string& failure_message) const;

    bool monitor_energy_;
    bool monitor_lz_;
    InvariantEvaluator carter_evaluator_;
    PhaseSpaceState initial_state_;
    double initial_energy_;
    double initial_lz_;
    double initial_carter_ = 0.0;
};

} // namespace solar::relativity::detail
