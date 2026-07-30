#pragma once

#include "kerr_separated_state.h"

#include "solar/relativity/kerr_separated.h"

namespace solar::relativity::detail {

struct KerrSeparatedStepAttempt {
    numerics::Dopri5StepResult<5> step;
    bool radial_forbidden;
    bool polar_forbidden;
    bool invalid_metric_point;
    bool non_finite_stage;
    double first_radial_forbidden_coordinate;
    double first_polar_forbidden_coordinate;
};

KerrSeparatedStepAttempt attempt_kerr_separated_step(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& current,
    double current_mino,
    double attempted_step,
    const KerrSeparatedConfig& config);

} // namespace solar::relativity::detail
