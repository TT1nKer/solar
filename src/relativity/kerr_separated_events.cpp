#include "kerr_separated_events.h"

#include <cmath>
#include <exception>
#include <limits>
#include <utility>

namespace solar::relativity::detail {
namespace {

struct EventEvaluation {
    KerrSeparatedState separated_state;
    PhaseSpaceState public_state;
    double value;
};

struct LocatedEvent {
    KerrSeparatedEventStatus status;
    std::optional<KerrSeparatedEventHit> hit;
    std::string message;
};

KerrSeparatedEventSelection no_selection() {
    return KerrSeparatedEventSelection{
        KerrSeparatedEventStatus::None,
        std::nullopt,
        TerminationReason::UserEvent,
        "no event found",
    };
}

KerrSeparatedEventSelection failed_selection(
    std::string message) {
    return KerrSeparatedEventSelection{
        KerrSeparatedEventStatus::Failed,
        std::nullopt,
        TerminationReason::EventRootFailure,
        std::move(message),
    };
}

bool brackets_zero(double left, double right) noexcept {
    return std::signbit(left) != std::signbit(right);
}

bool direction_matches(
    EventDirection direction,
    double start_value,
    double end_value) noexcept {
    if (direction == EventDirection::Any) {
        return true;
    }
    if (direction == EventDirection::Increasing) {
        return end_value > start_value;
    }
    return end_value < start_value;
}

template <std::size_t N, typename Projector>
EventEvaluation evaluate_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const numerics::Dopri5DenseOutput<N>& dense_output,
    const GeodesicEvent& event,
    double mino_parameter,
    const Projector& projector) {
    KerrSeparatedState state =
        projector(dense_output.evaluate(mino_parameter));
    PhaseSpaceState public_state =
        reconstruct_kerr_phase_space(
            metric, constants, state);
    const double value = event.function(public_state);
    if (!std::isfinite(value)) {
        throw std::domain_error(
            "event value is non-finite");
    }
    return EventEvaluation{
        std::move(state),
        std::move(public_state),
        value,
    };
}

KerrSeparatedEventHit make_hit(
    std::size_t event_index,
    double mino_parameter,
    EventEvaluation evaluation,
    std::size_t iterations) {
    EventHit public_hit{
        event_index,
        evaluation.public_state.affine,
        std::move(evaluation.public_state),
        evaluation.value,
        iterations,
    };
    return KerrSeparatedEventHit{
        std::move(public_hit),
        std::move(evaluation.separated_state),
        mino_parameter,
    };
}

template <std::size_t N, typename Projector>
LocatedEvent locate_event(
    std::size_t event_index,
    const GeodesicEvent& event,
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const numerics::Dopri5DenseOutput<N>& dense_output,
    const Projector& projector) {
    const double start_mino = dense_output.start();
    const double end_mino = dense_output.end();
    const double span = end_mino - start_mino;
    if (!std::isfinite(span) || span == 0.0) {
        return LocatedEvent{
            KerrSeparatedEventStatus::Failed,
            std::nullopt,
            "dense Mino interval is invalid",
        };
    }

    EventEvaluation start;
    EventEvaluation end;
    try {
        start = evaluate_event(
            metric,
            constants,
            dense_output,
            event,
            start_mino,
            projector);
        end = evaluate_event(
            metric,
            constants,
            dense_output,
            event,
            end_mino,
            projector);
    } catch (const std::exception& error) {
        return LocatedEvent{
            KerrSeparatedEventStatus::Failed,
            std::nullopt,
            std::string("event endpoint evaluation failed: ") +
                error.what(),
        };
    }

    if (!direction_matches(
            event.direction, start.value, end.value)) {
        return LocatedEvent{
            KerrSeparatedEventStatus::None,
            std::nullopt,
            "event direction does not match",
        };
    }
    if (start.value == 0.0) {
        return LocatedEvent{
            KerrSeparatedEventStatus::Found,
            make_hit(
                event_index,
                start_mino,
                std::move(start),
                0),
            "event root found",
        };
    }
    if (end.value == 0.0) {
        return LocatedEvent{
            KerrSeparatedEventStatus::Found,
            make_hit(
                event_index,
                end_mino,
                std::move(end),
                0),
            "event root found",
        };
    }
    if (!brackets_zero(start.value, end.value)) {
        return LocatedEvent{
            KerrSeparatedEventStatus::None,
            std::nullopt,
            "event does not bracket zero",
        };
    }

    double left = 0.0;
    double right = 1.0;
    EventEvaluation left_evaluation = std::move(start);
    EventEvaluation right_evaluation = std::move(end);

    for (std::size_t iteration = 1;
         iteration <= 100;
         ++iteration) {
        if (std::fabs(
                right_evaluation.public_state.affine -
                left_evaluation.public_state.affine) <=
            event.root_tolerance) {
            if (std::fabs(left_evaluation.value) <=
                std::fabs(right_evaluation.value)) {
                const double mino =
                    start_mino + span * left;
                return LocatedEvent{
                    KerrSeparatedEventStatus::Found,
                    make_hit(
                        event_index,
                        mino,
                        std::move(left_evaluation),
                        iteration - 1),
                    "event root found",
                };
            }
            const double mino =
                start_mino + span * right;
            return LocatedEvent{
                KerrSeparatedEventStatus::Found,
                make_hit(
                    event_index,
                    mino,
                    std::move(right_evaluation),
                    iteration - 1),
                "event root found",
            };
        }

        const double width = right - left;
        const double denominator =
            right_evaluation.value -
            left_evaluation.value;
        double candidate =
            right -
            right_evaluation.value * width / denominator;
        const double safeguard_left =
            left + 0.1 * width;
        const double safeguard_right =
            right - 0.1 * width;
        if (!std::isfinite(candidate) ||
            candidate <= safeguard_left ||
            candidate >= safeguard_right) {
            candidate = left + 0.5 * width;
        }

        const double mino =
            start_mino + span * candidate;
        EventEvaluation evaluation;
        try {
            evaluation = evaluate_event(
                metric,
                constants,
                dense_output,
                event,
                mino,
                projector);
        } catch (const std::exception& error) {
            return LocatedEvent{
                KerrSeparatedEventStatus::Failed,
                std::nullopt,
                std::string("event root evaluation failed: ") +
                    error.what(),
            };
        }
        if (evaluation.value == 0.0) {
            return LocatedEvent{
                KerrSeparatedEventStatus::Found,
                make_hit(
                    event_index,
                    mino,
                    std::move(evaluation),
                    iteration),
                "event root found",
            };
        }
        if (brackets_zero(
                left_evaluation.value,
                evaluation.value)) {
            right = candidate;
            right_evaluation = std::move(evaluation);
        } else {
            left = candidate;
            left_evaluation = std::move(evaluation);
        }
    }

    return LocatedEvent{
        KerrSeparatedEventStatus::Failed,
        std::nullopt,
        "event root iteration limit exceeded",
    };
}

template <std::size_t N, typename Projector>
KerrSeparatedEventSelection select_first_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const numerics::Dopri5DenseOutput<N>& dense_output,
    const std::vector<GeodesicEvent>& events,
    const Projector& projector) {
    std::optional<KerrSeparatedEventHit> selected;
    TerminationReason selected_reason =
        TerminationReason::UserEvent;
    double selected_fraction =
        std::numeric_limits<double>::infinity();
    const double span =
        dense_output.end() - dense_output.start();

    for (std::size_t index = 0;
         index < events.size();
         ++index) {
        const LocatedEvent located = locate_event(
            index,
            events[index],
            metric,
            constants,
            dense_output,
            projector);
        if (located.status ==
            KerrSeparatedEventStatus::Failed) {
            return failed_selection(located.message);
        }
        if (located.status !=
            KerrSeparatedEventStatus::Found) {
            continue;
        }

        const double fraction =
            (located.hit->mino_parameter -
             dense_output.start()) /
            span;
        if (fraction < selected_fraction) {
            selected_fraction = fraction;
            selected = located.hit;
            selected_reason = events[index].reason;
        }
    }

    if (!selected.has_value()) {
        return no_selection();
    }
    return KerrSeparatedEventSelection{
        KerrSeparatedEventStatus::Found,
        std::move(selected),
        selected_reason,
        "geodesic event reached",
    };
}

} // namespace

KerrSeparatedEventSelection select_first_kerr_step_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& directions,
    const numerics::Dopri5DenseOutput<5>& dense_output,
    const std::vector<GeodesicEvent>& events) {
    const auto projector =
        [&directions](const numerics::StateN<5>& values) {
            return KerrSeparatedState{
                values,
                directions.radial_direction,
                directions.polar_direction,
            };
        };
    return select_first_event(
        metric,
        constants,
        dense_output,
        events,
        projector);
}

KerrSeparatedEventSelection select_first_kerr_phase_step_event(
    const KerrBoyerLindquistMetric& metric,
    const KerrConstants& constants,
    const KerrSeparatedState& direction_fallback,
    const numerics::Dopri5DenseOutput<7>& dense_output,
    const std::vector<GeodesicEvent>& events) {
    const auto projector =
        [&direction_fallback](
            const KerrTurningPhaseState& phase) {
            return project_kerr_turning_phase(
                phase, direction_fallback);
        };
    return select_first_event(
        metric,
        constants,
        dense_output,
        events,
        projector);
}

} // namespace solar::relativity::detail
