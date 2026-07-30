#pragma once

#include "kerr_separated_turning.h"

#include "solar/relativity/geodesic_types.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace solar::relativity::detail {

enum class KerrPhaseFlowStatus {
    Accepted,
    Rejected,
    Terminated,
};

struct KerrPhaseFlowResult {
    KerrPhaseFlowStatus status;
    KerrTurningPhaseState phase_state;
    KerrSeparatedState separated_state;
    PhaseSpaceState public_state;
    double accepted_mino;
    double next_step;
    double additional_minimum_radius_M;
    std::size_t radial_turns;
    std::size_t polar_turns;
    bool non_finite_rejection;
    TerminationReason reason;
    std::string message;
    std::optional<EventHit> event;
};

KerrPhaseFlowResult advance_authoritative_kerr_phase(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedPotentials& potentials,
    const KerrSeparatedState& current,
    const PhaseSpaceState& current_public,
    const KerrTurningPhaseState& current_phase,
    double current_mino,
    double attempted_step,
    const KerrSeparatedConfig& config,
    const std::vector<GeodesicEvent>& events);

} // namespace solar::relativity::detail
