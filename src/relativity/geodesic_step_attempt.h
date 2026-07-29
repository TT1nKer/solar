#pragma once

#include "solar/relativity/geodesic_integrator.h"

namespace solar::relativity::detail {

struct GeodesicStepAttempt {
    numerics::Dopri5StepResult<8> step;
    bool invalid_metric_point;
    bool non_finite_stage;
};

GeodesicStepAttempt attempt_geodesic_step(
    const Metric& metric,
    const PhaseSpaceState& current,
    double attempted_step,
    const numerics::Dopri5Config<8>& config);

} // namespace solar::relativity::detail
