#pragma once

#include "solar/numerics/dopri5.h"
#include "solar/relativity/geodesic_types.h"

#include <optional>
#include <string>

namespace solar::relativity {

enum class EventRootStatus {
    NoRoot,
    Found,
    Failed,
};

struct EventRootResult {
    EventRootStatus status;
    std::optional<EventHit> hit;
    std::string message;
};

EventRootResult locate_event(
    std::size_t event_index,
    const GeodesicEvent& event,
    const numerics::Dopri5DenseOutput<8>& dense_output);

} // namespace solar::relativity
