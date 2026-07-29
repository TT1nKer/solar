#include "solar/relativity/geodesic_integrator.h"

#include "geodesic_config_internal.h"
#include "geodesic_event_selection.h"
#include "geodesic_step_attempt.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace solar::relativity {
namespace {

enum class RejectionKind {
    None,
    ErrorEstimate,
    InvalidMetricPoint,
    NonFiniteStage,
};

bool finite_phase_space(const PhaseSpaceState& state) {
    return std::isfinite(state.affine) &&
           state.x.v.all_finite() &&
           state.p.v.all_finite();
}

double invariant_drift(double current, double initial) {
    const double difference = std::fabs(current - initial);
    return initial == 0.0
               ? difference
               : difference / std::fabs(initial);
}

TerminationReason rejection_reason(RejectionKind kind) {
    if (kind == RejectionKind::InvalidMetricPoint) {
        return TerminationReason::InvalidMetricPoint;
    }
    if (kind == RejectionKind::NonFiniteStage) {
        return TerminationReason::NonFiniteState;
    }
    return TerminationReason::StepUnderflow;
}

std::string rejection_message(
    RejectionKind kind,
    const char* suffix) {
    if (kind == RejectionKind::InvalidMetricPoint) {
        return std::string("metric-domain trial rejection ") + suffix;
    }
    if (kind == RejectionKind::NonFiniteStage) {
        return std::string("non-finite trial rejection ") + suffix;
    }
    return std::string("error-estimate rejection ") + suffix;
}

} // namespace

GeodesicIntegrator::GeodesicIntegrator(
    const Metric& metric) noexcept
    : metric_(&metric) {}

GeodesicIntegrationResult GeodesicIntegrator::integrate(
    const PhaseSpaceState& initial,
    const GeodesicIntegrationConfig& config,
    const std::vector<GeodesicEvent>& events) const {
    detail::validate_geodesic_config(config);

    IntegrationDiagnostics diagnostics{};
    const auto terminate =
        [&diagnostics](
            const PhaseSpaceState& state,
            TerminationReason reason,
            std::string message,
            std::optional<EventHit> event = std::nullopt) {
            diagnostics.reason = reason;
            diagnostics.message = std::move(message);
            return GeodesicIntegrationResult{
                state, diagnostics, std::move(event)};
        };

    if (!finite_phase_space(initial)) {
        return terminate(
            initial,
            TerminationReason::NonFiniteState,
            "initial phase-space state is non-finite");
    }
    if (!metric_->valid_point(initial.x)) {
        return terminate(
            initial,
            TerminationReason::InvalidMetricPoint,
            "initial state is outside the metric domain");
    }

    double initial_constraint;
    try {
        initial_constraint = hamiltonian_constraint_error(
            *metric_, initial, config.kind);
    } catch (const std::domain_error& error) {
        return terminate(
            initial,
            TerminationReason::NonFiniteState,
            std::string("initial Hamiltonian evaluation failed: ") +
                error.what());
    }
    diagnostics.max_constraint_error = initial_constraint;
    if (initial_constraint > config.constraint_tolerance) {
        return terminate(
            initial,
            TerminationReason::ConstraintViolation,
            "initial Hamiltonian constraint exceeds tolerance");
    }

    const double initial_energy = -initial.p.v[0];
    const double initial_lz = initial.p.v[3];
    if (config.monitor_energy) {
        diagnostics.max_energy_rel_error = 0.0;
    }
    if (config.monitor_lz) {
        diagnostics.max_lz_rel_error = 0.0;
    }

    std::vector<GeodesicEvent> active_events = events;
    if (std::isfinite(config.max_proper_time)) {
        active_events.push_back(GeodesicEvent{
            "proper-time limit",
            [origin = initial.affine,
             limit = config.max_proper_time](
                const PhaseSpaceState& state) {
                return std::fabs(state.affine - origin) - limit;
            },
            EventDirection::Increasing,
            TerminationReason::MaxProperTime,
            config.min_step,
        });
    }
    if (std::isfinite(config.max_coordinate_time)) {
        active_events.push_back(GeodesicEvent{
            "coordinate-time limit",
            [origin = initial.x.v[0],
             limit = config.max_coordinate_time](
                const PhaseSpaceState& state) {
                return std::fabs(state.x.v[0] - origin) - limit;
            },
            EventDirection::Any,
            TerminationReason::MaxCoordinateTime,
            config.min_step,
        });
    }

    const detail::GeodesicEventSelection initial_event =
        detail::select_initial_any_event(initial, active_events);
    if (initial_event.status ==
        detail::GeodesicEventSelectionStatus::Failed) {
        return terminate(
            initial,
            TerminationReason::EventRootFailure,
            initial_event.message);
    }
    if (initial_event.status ==
        detail::GeodesicEventSelectionStatus::Found) {
        return terminate(
            initial,
            initial_event.reason,
            initial_event.message,
            initial_event.hit);
    }

    PhaseSpaceState current = initial;
    double proposed_step = config.initial_step;
    const double direction =
        std::copysign(1.0, config.initial_step);
    std::size_t attempted_steps = 0;
    std::size_t consecutive_rejections = 0;
    RejectionKind last_rejection = RejectionKind::None;

    while (true) {
        if (attempted_steps >= config.max_total_steps) {
            return terminate(
                current,
                TerminationReason::MaxSteps,
                "maximum total DOPRI5 attempts reached");
        }

        const double affine_displacement =
            std::fabs(current.affine - initial.affine);
        const double affine_remaining =
            config.max_affine - affine_displacement;
        if (affine_remaining <= 0.0) {
            return terminate(
                current,
                TerminationReason::MaxAffine,
                "maximum affine displacement reached");
        }

        const double requested_magnitude = std::min(
            {std::fabs(proposed_step),
             config.max_step,
             affine_remaining});
        const bool lands_on_affine_limit =
            requested_magnitude == affine_remaining;
        if (!std::isfinite(requested_magnitude) ||
            requested_magnitude == 0.0 ||
            (requested_magnitude < config.min_step &&
             !lands_on_affine_limit)) {
            return terminate(
                current,
                rejection_reason(last_rejection),
                last_rejection == RejectionKind::None
                    ? "proposed step is below the configured minimum"
                    : rejection_message(
                          last_rejection,
                          "fell below the configured minimum"));
        }
        const double attempted_step =
            direction * requested_magnitude;

        ++attempted_steps;
        const detail::GeodesicStepAttempt trial =
            detail::attempt_geodesic_step(
                *metric_,
                current,
                attempted_step,
                config.dopri5);
        const auto& step = trial.step;

        if (step.status !=
            numerics::Dopri5StepResult<8>::Status::Completed ||
            !step.accepted) {
            ++diagnostics.rejected_steps;
            ++consecutive_rejections;
            if (trial.invalid_metric_point) {
                last_rejection =
                    RejectionKind::InvalidMetricPoint;
                proposed_step =
                    direction * requested_magnitude *
                    config.dopri5.min_factor;
            } else if (
                trial.non_finite_stage ||
                step.status !=
                    numerics::Dopri5StepResult<8>::Status::Completed) {
                last_rejection = RejectionKind::NonFiniteStage;
                proposed_step =
                    direction * requested_magnitude *
                    config.dopri5.min_factor;
            } else {
                last_rejection = RejectionKind::ErrorEstimate;
                proposed_step = step.next_step;
            }

            if (consecutive_rejections >=
                config.max_rejections_per_step) {
                return terminate(
                    current,
                    rejection_reason(last_rejection),
                    rejection_message(
                        last_rejection,
                        "reached the rejection limit"));
            }
            continue;
        }

        if (!step.dense_output.has_value()) {
            return terminate(
                current,
                TerminationReason::NonFiniteState,
                "accepted DOPRI5 step has no dense output");
        }

        const detail::GeodesicEventSelection selected_event =
            detail::select_first_step_event(
                *step.dense_output, active_events);
        if (selected_event.status ==
            detail::GeodesicEventSelectionStatus::Failed) {
            return terminate(
                current,
                TerminationReason::EventRootFailure,
                selected_event.message);
        }

        PhaseSpaceState accepted_state;
        double accepted_step_magnitude;
        if (selected_event.hit.has_value()) {
            accepted_state = selected_event.hit->state;
            accepted_step_magnitude = std::fabs(
                accepted_state.affine - current.affine);
        } else {
            accepted_state = unpack_phase_space(
                current.affine + attempted_step,
                step.state);
            accepted_step_magnitude =
                std::fabs(attempted_step);
        }

        if (selected_event.hit.has_value() &&
            accepted_step_magnitude == 0.0) {
            return terminate(
                accepted_state,
                selected_event.reason,
                "geodesic event is at the current state",
                selected_event.hit);
        }

        if (!finite_phase_space(accepted_state)) {
            return terminate(
                current,
                TerminationReason::NonFiniteState,
                "accepted phase-space state is non-finite");
        }
        if (!metric_->valid_point(accepted_state.x)) {
            return terminate(
                current,
                TerminationReason::InvalidMetricPoint,
                "accepted state is outside the metric domain");
        }

        double constraint;
        try {
            constraint = hamiltonian_constraint_error(
                *metric_, accepted_state, config.kind);
        } catch (const std::domain_error& error) {
            return terminate(
                current,
                TerminationReason::NonFiniteState,
                std::string(
                    "accepted Hamiltonian evaluation failed: ") +
                    error.what());
        }

        ++diagnostics.accepted_steps;
        diagnostics.min_step =
            std::isnan(diagnostics.min_step)
                ? accepted_step_magnitude
                : std::min(
                      diagnostics.min_step,
                      accepted_step_magnitude);
        diagnostics.max_step =
            std::isnan(diagnostics.max_step)
                ? accepted_step_magnitude
                : std::max(
                      diagnostics.max_step,
                      accepted_step_magnitude);
        diagnostics.max_constraint_error = std::max(
            diagnostics.max_constraint_error, constraint);
        if (config.monitor_energy) {
            diagnostics.max_energy_rel_error = std::max(
                diagnostics.max_energy_rel_error,
                invariant_drift(
                    -accepted_state.p.v[0],
                    initial_energy));
        }
        if (config.monitor_lz) {
            diagnostics.max_lz_rel_error = std::max(
                diagnostics.max_lz_rel_error,
                invariant_drift(
                    accepted_state.p.v[3],
                    initial_lz));
        }
        current = accepted_state;

        if (constraint > config.constraint_tolerance) {
            return terminate(
                current,
                TerminationReason::ConstraintViolation,
                "accepted Hamiltonian constraint exceeds tolerance");
        }
        if (selected_event.hit.has_value()) {
            return terminate(
                current,
                selected_event.reason,
                selected_event.message,
                selected_event.hit);
        }

        consecutive_rejections = 0;
        last_rejection = RejectionKind::None;
        proposed_step = step.next_step;

        if (std::fabs(current.affine - initial.affine) >=
            config.max_affine) {
            return terminate(
                current,
                TerminationReason::MaxAffine,
                "maximum affine displacement reached");
        }
    }
}

} // namespace solar::relativity
