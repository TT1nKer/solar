#include "kerr_separated_phase_flow.h"

#include "kerr_separated_events.h"
#include "kerr_separated_state.h"

#include <cmath>
#include <exception>
#include <initializer_list>
#include <limits>
#include <string>
#include <utility>

namespace solar::relativity::detail {
namespace {

bool finite_phase_space(const PhaseSpaceState& state) noexcept {
    return std::isfinite(state.affine) &&
           state.x.v.all_finite() &&
           state.p.v.all_finite();
}

KerrPhaseFlowResult terminated(
    const KerrTurningPhaseState& phase,
    const KerrSeparatedState& separated,
    PhaseSpaceState public_state,
    double current_mino,
    TerminationReason reason,
    std::string message) {
    return KerrPhaseFlowResult{
        KerrPhaseFlowStatus::Terminated,
        phase,
        separated,
        std::move(public_state),
        current_mino,
        0.0,
        std::numeric_limits<double>::quiet_NaN(),
        0,
        0,
        false,
        reason,
        std::move(message),
        std::nullopt,
    };
}

bool coordinate_is_locked(
    TurningCoordinate coordinate,
    const KerrSeparatedState& state) noexcept {
    return coordinate == TurningCoordinate::Radial
               ? state.radial_direction ==
                     SeparatedDirection::Locked
               : state.polar_direction ==
                     SeparatedDirection::Locked;
}

} // namespace

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
    const std::vector<GeodesicEvent>& events) {
    const auto step = attempt_kerr_turning_phase_step(
        metric,
        constants,
        potentials,
        current_phase,
        current_mino,
        attempted_step,
        config);
    if (step.status !=
            numerics::Dopri5StepResult<7>::Status::Completed ||
        !step.accepted) {
        return KerrPhaseFlowResult{
            KerrPhaseFlowStatus::Rejected,
            current_phase,
            current,
            current_public,
            current_mino,
            step.status ==
                    numerics::Dopri5StepResult<7>::
                        Status::Completed
                ? step.next_step
                : attempted_step *
                      config.dopri5.min_factor,
            std::numeric_limits<double>::quiet_NaN(),
            0,
            0,
            step.status !=
                numerics::Dopri5StepResult<7>::
                    Status::Completed,
            TerminationReason::StepUnderflow,
            "turning phase step rejected",
            std::nullopt,
        };
    }
    if (!step.dense_output.has_value()) {
        return terminated(
            current_phase,
            current,
            current_public,
            current_mino,
            TerminationReason::NonFiniteState,
            "accepted turning phase step has no dense output");
    }

    const auto selected_event =
        select_first_kerr_phase_step_event(
            metric,
            constants,
            current,
            *step.dense_output,
            events);
    if (selected_event.status ==
        KerrSeparatedEventStatus::Failed) {
        return terminated(
            current_phase,
            current,
            current_public,
            current_mino,
            TerminationReason::EventRootFailure,
            selected_event.message);
    }

    const double accepted_mino =
        selected_event.hit.has_value()
            ? selected_event.hit->mino_parameter
            : current_mino + attempted_step;
    const auto accepted_phase =
        step.dense_output->evaluate(accepted_mino);
    const KerrSeparatedState accepted =
        project_kerr_turning_phase(
            accepted_phase, current);
    PhaseSpaceState accepted_public;
    std::optional<EventHit> public_hit;
    if (selected_event.hit.has_value()) {
        accepted_public =
            selected_event.hit->public_hit.state;
        public_hit = selected_event.hit->public_hit;
    } else {
        try {
            accepted_public = reconstruct_kerr_phase_space(
                metric, constants, accepted);
        } catch (const std::exception& error) {
            return terminated(
                current_phase,
                current,
                current_public,
                current_mino,
                TerminationReason::NonFiniteState,
                std::string(
                    "accepted turning phase state failed: ") +
                    error.what());
        }
    }
    if (!finite_phase_space(accepted_public) ||
        !metric.valid_point(accepted_public.x)) {
        return terminated(
            current_phase,
            current,
            current_public,
            current_mino,
            TerminationReason::InvalidMetricPoint,
            "accepted turning phase state is outside Kerr BL");
    }

    double additional_minimum_radius =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t radial_turns = 0;
    std::size_t polar_turns = 0;
    for (const TurningCoordinate coordinate :
         {TurningCoordinate::Radial,
          TurningCoordinate::Polar}) {
        if (coordinate_is_locked(coordinate, current)) {
            continue;
        }
        const auto crossing =
            locate_kerr_turning_phase_crossing(
                coordinate,
                potentials,
                *step.dense_output,
                accepted_mino,
                config.potential_tolerance,
                config.critical_derivative_tolerance);
        if (!crossing.has_value()) {
            continue;
        }
        if (crossing->status == TurningStatus::Failed) {
            return terminated(
                current_phase,
                current,
                current_public,
                current_mino,
                TerminationReason::EventRootFailure,
                crossing->message);
        }
        if (crossing->status ==
            TurningStatus::NearCritical) {
            const auto critical_phase =
                step.dense_output->evaluate(
                    crossing->mino_parameter);
            const auto critical_separated =
                project_kerr_turning_phase(
                    critical_phase, current);
            PhaseSpaceState critical_public = current_public;
            try {
                critical_public =
                    reconstruct_kerr_phase_space(
                        metric,
                        constants,
                        critical_separated);
            } catch (const std::exception&) {
                // The last valid public state remains the safe result.
            }
            return terminated(
                critical_phase,
                critical_separated,
                std::move(critical_public),
                crossing->mino_parameter,
                TerminationReason::NearCriticalOrbit,
                crossing->message);
        }
        if (coordinate == TurningCoordinate::Radial) {
            ++radial_turns;
            additional_minimum_radius =
                crossing->root_radius_M;
        } else {
            ++polar_turns;
        }
    }

    return KerrPhaseFlowResult{
        KerrPhaseFlowStatus::Accepted,
        accepted_phase,
        accepted,
        std::move(accepted_public),
        accepted_mino,
        step.next_step,
        additional_minimum_radius,
        radial_turns,
        polar_turns,
        false,
        selected_event.reason,
        selected_event.message,
        std::move(public_hit),
    };
}

} // namespace solar::relativity::detail
