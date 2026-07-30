#pragma once

#include "kerr_separated_state.h"

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

} // namespace solar::relativity::detail
