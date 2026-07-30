#include "geodesic_event_selection.h"

#include <cmath>
#include <exception>
#include <limits>
#include <utility>

namespace solar::relativity::detail {
namespace {

GeodesicEventSelection no_event() {
    return GeodesicEventSelection{
        GeodesicEventSelectionStatus::None,
        std::nullopt,
        TerminationReason::UserEvent,
        "no event found",
    };
}

GeodesicEventSelection failed_event(std::string message) {
    return GeodesicEventSelection{
        GeodesicEventSelectionStatus::Failed,
        std::nullopt,
        TerminationReason::EventRootFailure,
        std::move(message),
    };
}

} // namespace

std::optional<std::string> validate_geodesic_event_contracts(
    const std::vector<GeodesicEvent>& events) {
    for (const GeodesicEvent& event : events) {
        if (!event.function ||
            !is_valid_event_direction(event.direction) ||
            !std::isfinite(event.root_tolerance) ||
            event.root_tolerance <= 0.0) {
            return "initial event contract is invalid";
        }
    }
    return std::nullopt;
}

GeodesicEventSelection select_initial_any_event(
    const PhaseSpaceState& initial,
    const std::vector<GeodesicEvent>& events) {
    for (std::size_t index = 0; index < events.size(); ++index) {
        const GeodesicEvent& event = events[index];
        if (event.direction != EventDirection::Any) {
            continue;
        }

        double value;
        try {
            value = event.function(initial);
        } catch (const std::exception& error) {
            return failed_event(
                std::string("initial event evaluation failed: ") +
                error.what());
        }
        if (!std::isfinite(value)) {
            return failed_event(
                "initial event value is non-finite");
        }
        if (value == 0.0) {
            return GeodesicEventSelection{
                GeodesicEventSelectionStatus::Found,
                EventHit{
                    index,
                    initial.affine,
                    initial,
                    value,
                    0,
                },
                event.reason,
                "geodesic event is at the initial state",
            };
        }
    }
    return no_event();
}

GeodesicEventSelection select_first_step_event(
    const numerics::Dopri5DenseOutput<8>& dense_output,
    const std::vector<GeodesicEvent>& events) {
    std::optional<EventHit> selected;
    TerminationReason selected_reason =
        TerminationReason::UserEvent;
    double selected_fraction =
        std::numeric_limits<double>::infinity();
    const double span = dense_output.end() - dense_output.start();

    for (std::size_t index = 0; index < events.size(); ++index) {
        const EventRootResult root =
            locate_event(index, events[index], dense_output);
        if (root.status == EventRootStatus::Failed) {
            return failed_event(root.message);
        }
        if (root.status != EventRootStatus::Found) {
            continue;
        }
        const double fraction =
            (root.hit->affine - dense_output.start()) / span;
        if (fraction < selected_fraction) {
            selected_fraction = fraction;
            selected = root.hit;
            selected_reason = events[index].reason;
        }
    }

    if (!selected.has_value()) {
        return no_event();
    }
    return GeodesicEventSelection{
        GeodesicEventSelectionStatus::Found,
        selected,
        selected_reason,
        "geodesic event reached",
    };
}

} // namespace solar::relativity::detail
