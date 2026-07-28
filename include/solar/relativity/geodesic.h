#pragma once

#include "solar/relativity/schwarzschild.h"
#include <vector>

namespace solar::relativity {

enum class PropagationStatus {
    Completed,
    HorizonReached,
    StepLimitReached,
    StepSizeUnderflow,
};

struct GeodesicSample {
    double affine_parameter = 0.0;
    GeodesicState state;
    double metric_norm = 0.0;
    ConservedQuantities conserved;
};

struct GeodesicOptions {
    double affine_duration = 1.0;
    double initial_step = 1e-3;
    double absolute_tolerance = 1e-11;
    double relative_tolerance = 1e-11;
    double horizon_margin = 1e-6;
    int max_steps = 100000;
};

struct GeodesicResult {
    PropagationStatus status = PropagationStatus::Completed;
    std::vector<GeodesicSample> samples;
    int accepted_steps = 0;
    int rejected_steps = 0;
    double max_metric_norm_drift = 0.0;
    double max_relative_energy_drift = 0.0;
    double max_relative_angular_momentum_drift = 0.0;
};

GeodesicResult propagate_geodesic(
    const SchwarzschildSpacetime& spacetime,
    const GeodesicState& initial_state,
    GeodesicKind kind,
    const GeodesicOptions& options);

} // namespace solar::relativity
