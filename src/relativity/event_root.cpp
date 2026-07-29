#include "solar/relativity/event_root.h"

#include "solar/relativity/hamiltonian.h"

#include <algorithm>
#include <cmath>
#include <exception>

namespace solar::relativity {
namespace {

EventRootResult no_root(const std::string& message) {
    return EventRootResult{
        EventRootStatus::NoRoot, std::nullopt, message};
}

EventRootResult failed_root(const std::string& message) {
    return EventRootResult{
        EventRootStatus::Failed, std::nullopt, message};
}

EventRootResult found_root(
    std::size_t event_index,
    double affine,
    const numerics::StateN<8>& packed,
    double value,
    std::size_t iterations) {
    return EventRootResult{
        EventRootStatus::Found,
        EventHit{
            event_index,
            affine,
            unpack_phase_space(affine, packed),
            value,
            iterations,
        },
        "event root found",
    };
}

bool direction_matches(
    EventDirection direction,
    double start_value,
    double end_value) {
    if (direction == EventDirection::Any) {
        return true;
    }
    if (direction == EventDirection::Increasing) {
        return end_value > start_value;
    }
    return end_value < start_value;
}

bool brackets_zero(double left, double right) {
    return std::signbit(left) != std::signbit(right);
}

} // namespace

EventRootResult locate_event(
    std::size_t event_index,
    const GeodesicEvent& event,
    const numerics::Dopri5DenseOutput<8>& dense_output) {
    if (!event.function) {
        return failed_root("event function is empty");
    }
    if (!std::isfinite(event.root_tolerance) ||
        event.root_tolerance <= 0.0) {
        return failed_root(
            "event root tolerance must be positive and finite");
    }

    const double affine_start = dense_output.start();
    const double affine_end = dense_output.end();
    const double affine_span = affine_end - affine_start;
    if (!std::isfinite(affine_span) || affine_span == 0.0) {
        return failed_root(
            "dense-output affine interval is invalid");
    }

    numerics::StateN<8> start_packed;
    numerics::StateN<8> end_packed;
    double start_value;
    double end_value;
    try {
        start_packed = dense_output.evaluate(affine_start);
        end_packed = dense_output.evaluate(affine_end);
        start_value = event.function(
            unpack_phase_space(affine_start, start_packed));
        end_value = event.function(
            unpack_phase_space(affine_end, end_packed));
    } catch (const std::exception& error) {
        return failed_root(
            std::string("event endpoint evaluation failed: ") +
            error.what());
    }
    if (!std::isfinite(start_value) ||
        !std::isfinite(end_value)) {
        return failed_root(
            "event endpoint value is non-finite");
    }
    if (!direction_matches(
            event.direction, start_value, end_value)) {
        return no_root("event direction does not match");
    }
    if (start_value == 0.0) {
        return found_root(
            event_index,
            affine_start,
            start_packed,
            start_value,
            0);
    }
    if (end_value == 0.0) {
        return found_root(
            event_index,
            affine_end,
            end_packed,
            end_value,
            0);
    }
    if (!brackets_zero(start_value, end_value)) {
        return no_root("event values do not bracket zero");
    }

    double left = 0.0;
    double right = 1.0;
    double left_value = start_value;
    double right_value = end_value;
    numerics::StateN<8> left_packed = start_packed;
    numerics::StateN<8> right_packed = end_packed;

    for (std::size_t iteration = 1; iteration <= 100; ++iteration) {
        const double bracket_width = right - left;
        if (std::fabs(affine_span) * bracket_width <=
            event.root_tolerance) {
            if (std::fabs(left_value) <= std::fabs(right_value)) {
                const double affine =
                    affine_start + affine_span * left;
                return found_root(
                    event_index,
                    affine,
                    left_packed,
                    left_value,
                    iteration - 1);
            }
            const double affine =
                affine_start + affine_span * right;
            return found_root(
                event_index,
                affine,
                right_packed,
                right_value,
                iteration - 1);
        }

        const double denominator = right_value - left_value;
        double candidate =
            right -
            right_value * bracket_width / denominator;
        const double safeguard_left =
            left + 0.1 * bracket_width;
        const double safeguard_right =
            right - 0.1 * bracket_width;
        if (!std::isfinite(candidate) ||
            candidate <= safeguard_left ||
            candidate >= safeguard_right) {
            candidate = left + 0.5 * bracket_width;
        }

        const double affine =
            affine_start + affine_span * candidate;
        numerics::StateN<8> packed;
        double value;
        try {
            packed = dense_output.evaluate(affine);
            value = event.function(
                unpack_phase_space(affine, packed));
        } catch (const std::exception& error) {
            return failed_root(
                std::string("event root evaluation failed: ") +
                error.what());
        }
        if (!std::isfinite(value)) {
            return failed_root(
                "event root value is non-finite");
        }
        if (value == 0.0) {
            return found_root(
                event_index,
                affine,
                packed,
                value,
                iteration);
        }

        if (brackets_zero(left_value, value)) {
            right = candidate;
            right_value = value;
            right_packed = packed;
        } else {
            left = candidate;
            left_value = value;
            left_packed = packed;
        }
    }

    return failed_root(
        "event root iteration limit exceeded");
}

} // namespace solar::relativity
