#pragma once

#include "solar/relativity/event_root.h"

#include <optional>
#include <string>
#include <vector>

namespace solar::relativity::detail {

enum class GeodesicEventSelectionStatus {
    None,
    Found,
    Failed,
};

struct GeodesicEventSelection {
    GeodesicEventSelectionStatus status;
    std::optional<EventHit> hit;
    TerminationReason reason;
    std::string message;
};

GeodesicEventSelection select_initial_any_event(
    const PhaseSpaceState& initial,
    const std::vector<GeodesicEvent>& events);

GeodesicEventSelection select_first_step_event(
    const numerics::Dopri5DenseOutput<8>& dense_output,
    const std::vector<GeodesicEvent>& events);

} // namespace solar::relativity::detail
