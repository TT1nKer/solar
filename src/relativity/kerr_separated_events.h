#pragma once

#include "kerr_separated_state.h"
#include "kerr_separated_turning.h"

#include "solar/relativity/geodesic_types.h"

#include <optional>
#include <string>
#include <vector>

namespace solar::relativity::detail {

enum class KerrSeparatedEventStatus {
    None,
    Found,
    Failed,
};

struct KerrSeparatedEventHit {
    EventHit public_hit;
    KerrSeparatedState separated_state;
    double mino_parameter;
};

struct KerrSeparatedEventSelection {
    KerrSeparatedEventStatus status;
    std::optional<KerrSeparatedEventHit> hit;
    TerminationReason reason;
    std::string message;
};

KerrSeparatedEventSelection select_first_kerr_step_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& directions,
    const numerics::Dopri5DenseOutput<5>& dense_output,
    const std::vector<GeodesicEvent>& events);

KerrSeparatedEventSelection select_first_kerr_phase_step_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& direction_fallback,
    const numerics::Dopri5DenseOutput<7>& dense_output,
    const std::vector<GeodesicEvent>& events);

KerrSeparatedEventSelection
select_first_kerr_phase_interval_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrSeparatedState& direction_fallback,
    const KerrTurningPhaseState& initial_phase,
    double initial_mino,
    double mino_step,
    const KerrSeparatedConfig& config,
    const std::vector<GeodesicEvent>& events);

KerrSeparatedEventSelection
select_first_kerr_turning_release_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrSeparatedState& direction_fallback,
    const KerrTurningPhaseState& initial_phase,
    double initial_mino,
    double release_mino_step,
    double root_mino_step,
    const KerrSeparatedConfig& config,
    const std::vector<GeodesicEvent>& events);

} // namespace solar::relativity::detail
