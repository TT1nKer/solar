#include "solar/relativity/kerr_separated.h"

#include "geodesic_event_selection.h"
#include "kerr_separated_config_internal.h"
#include "kerr_separated_diagnostics.h"
#include "kerr_separated_events.h"
#include "kerr_separated_state.h"
#include "kerr_separated_step.h"
#include "kerr_separated_turning.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

namespace solar::relativity {
namespace {

enum class RejectionKind {
    None,
    ErrorEstimate,
    ForbiddenPotential,
    InvalidMetricPoint,
    NonFiniteStage,
};

struct PendingTurning {
    detail::TurningCoordinate coordinate;
    detail::TurningRoot root;
};

bool finite_phase_space(const PhaseSpaceState& state) noexcept {
    return std::isfinite(state.affine) &&
           state.x.v.all_finite() &&
           state.p.v.all_finite();
}

TerminationReason rejection_reason(
    RejectionKind kind) noexcept {
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
    if (kind == RejectionKind::ForbiddenPotential) {
        return std::string("turning-potential trial rejection ") +
               suffix;
    }
    if (kind == RejectionKind::InvalidMetricPoint) {
        return std::string("metric-domain trial rejection ") +
               suffix;
    }
    if (kind == RejectionKind::NonFiniteStage) {
        return std::string("non-finite trial rejection ") +
               suffix;
    }
    return std::string("error-estimate rejection ") + suffix;
}

double affine_root_tolerance(
    const KerrBoyerLindquistMetric& metric,
    const KerrSeparatedConfig& config) noexcept {
    return config.root_tolerance * metric.mass();
}

bool ready_for_turning_release(
    const PendingTurning& pending,
    const detail::KerrSeparatedState& current,
    const detail::KerrSeparatedPotentials& potentials,
    const KerrSeparatedConfig& config) {
    const auto values = potentials.evaluate(
        current.values[detail::kRadius],
        current.values[detail::kMu]);
    const bool radial =
        pending.coordinate ==
        detail::TurningCoordinate::Radial;
    const double potential =
        radial ? values.radial : values.polar;
    const double potential_scale =
        radial ? values.radial_scale : values.polar_scale;
    const double coordinate =
        current.values[
            radial ? detail::kRadius : detail::kMu];
    const double coordinate_scale =
        radial ? potentials.mass() : 1.0;
    return std::fabs(potential) / potential_scale <=
               config.potential_tolerance ||
           std::fabs(
               coordinate - pending.root.coordinate) /
                   coordinate_scale <=
               config.root_tolerance;
}

} // namespace

KerrSeparatedIntegrator::KerrSeparatedIntegrator(
    const KerrBoyerLindquistMetric& metric) noexcept
    : metric_(&metric) {}

KerrSeparatedIntegrationResult
KerrSeparatedIntegrator::integrate(
    const PhaseSpaceState& initial,
    const KerrSeparatedConfig& config,
    const std::vector<GeodesicEvent>& events) const {
    detail::validate_kerr_separated_config(config);

    KerrSeparatedDiagnostics diagnostics{};
    KerrConstants constants{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
    };
    const auto terminate =
        [&diagnostics, &constants](
            const PhaseSpaceState& state,
            TerminationReason reason,
            std::string message,
            std::optional<EventHit> event = std::nullopt) {
            diagnostics.reason = reason;
            diagnostics.message = std::move(message);
            return KerrSeparatedIntegrationResult{
                state,
                constants,
                diagnostics,
                std::move(event),
            };
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
            "initial state is outside the Kerr BL domain");
    }

    std::vector<GeodesicEvent> active_events = events;
    active_events.push_back(GeodesicEvent{
        "affine-displacement limit",
        [origin = initial.affine,
         limit = config.max_affine](
            const PhaseSpaceState& state) {
            return std::fabs(state.affine - origin) - limit;
        },
        EventDirection::Any,
        TerminationReason::MaxAffine,
        affine_root_tolerance(*metric_, config),
    });
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
            affine_root_tolerance(*metric_, config),
        });
    }

    const std::optional<std::string> event_contract_failure =
        detail::validate_geodesic_event_contracts(active_events);
    if (event_contract_failure.has_value()) {
        return terminate(
            initial,
            TerminationReason::EventRootFailure,
            *event_contract_failure);
    }

    try {
        constants = evaluate_kerr_constants(
            *metric_, initial, config.kind);
    } catch (const std::exception& error) {
        return terminate(
            initial,
            TerminationReason::NonFiniteState,
            std::string(
                "initial Kerr constants evaluation failed: ") +
                error.what());
    }

    detail::KerrSeparatedInitialState initialized;
    try {
        initialized =
            detail::initialize_kerr_separated_state(
                *metric_,
                initial,
                config.kind,
                config.potential_tolerance,
                config.critical_derivative_tolerance,
                config.initial_mino_step);
    } catch (
        const detail::KerrSeparatedCriticalInitialState&
            error) {
        return terminate(
            initial,
            TerminationReason::NearCriticalOrbit,
            error.what());
    } catch (const std::domain_error& error) {
        return terminate(
            initial,
            TerminationReason::ConstraintViolation,
            error.what());
    }
    constants = initialized.constants;

    const double initial_mu =
        initialized.state.values[detail::kMu];
    if (1.0 - initial_mu * initial_mu <=
            config.polar_axis_tolerance &&
        constants.Lz != 0.0) {
        return terminate(
            initial,
            TerminationReason::InvalidMetricPoint,
            "near-axis nonzero-Lz Kerr BL motion is unsupported");
    }

    const detail::KerrSeparatedDiagnosticTracker
        diagnostic_tracker(
            *metric_, config.kind, constants, initial);
    try {
        diagnostic_tracker.initialize(
            initialized.state, initial, diagnostics);
    } catch (const std::exception& error) {
        return terminate(
            initial,
            TerminationReason::NonFiniteState,
            std::string(
                "initial separated diagnostics failed: ") +
                error.what());
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

    detail::KerrSeparatedState current =
        initialized.state;
    PhaseSpaceState current_public = initial;
    double current_mino = 0.0;
    double proposed_step = config.initial_mino_step;
    const double integration_direction =
        std::copysign(1.0, proposed_step);
    std::size_t attempted_steps = 0;
    std::size_t consecutive_rejections = 0;
    RejectionKind last_rejection = RejectionKind::None;
    std::optional<PendingTurning> pending_turning;
    const detail::KerrSeparatedPotentials potentials(
        metric_->mass(), metric_->spin_length(), constants);
    detail::KerrTurningPhaseState turning_phase;
    bool turning_phase_authoritative = false;
    try {
        turning_phase =
            detail::initialize_kerr_turning_phase(
                current, potentials);
    } catch (const std::exception& error) {
        return terminate(
            current_public,
            TerminationReason::NonFiniteState,
            std::string(
                "turning phase initialization failed: ") +
                error.what());
    }

    while (true) {
        if (attempted_steps >= config.max_total_steps) {
            return terminate(
                current_public,
                TerminationReason::MaxSteps,
                "maximum total separated DOPRI5 attempts reached");
        }

        if (pending_turning.has_value() &&
            ready_for_turning_release(
                *pending_turning,
                current,
                potentials,
                config)) {
            ++attempted_steps;
            detail::TurningRelease release;
            PhaseSpaceState released_public;
            try {
                release =
                    detail::release_kerr_turning_point(
                        pending_turning->coordinate,
                        pending_turning->root,
                        *metric_,
                        constants,
                        current,
                        turning_phase,
                        integration_direction,
                        config);
                released_public =
                    detail::reconstruct_kerr_phase_space(
                        *metric_, constants, release.state);
            } catch (const std::exception& error) {
                return terminate(
                    current_public,
                    TerminationReason::StepUnderflow,
                    std::string(
                        "turning release failed: ") +
                        error.what());
            }

            const double release_magnitude =
                std::fabs(release.mino_step);
            try {
                diagnostic_tracker.accept(
                    release.state,
                    released_public,
                    release_magnitude,
                    diagnostics,
                    release.root_radius_M);
            } catch (const std::exception& error) {
                return terminate(
                    current_public,
                    TerminationReason::NonFiniteState,
                    std::string(
                        "turning diagnostics failed: ") +
                        error.what());
            }
            if (pending_turning->coordinate ==
                detail::TurningCoordinate::Radial) {
                ++diagnostics.radial_turns;
            } else {
                ++diagnostics.polar_turns;
            }

            current = release.state;
            turning_phase = release.phase_state;
            turning_phase_authoritative = true;
            current_public = released_public;
            current_mino += release.mino_step;
            proposed_step = std::copysign(
                std::min(
                    config.max_mino_step,
                    std::max(
                        config.min_mino_step,
                        2.0 * release_magnitude)),
                integration_direction);
            consecutive_rejections = 0;
            last_rejection = RejectionKind::None;
            pending_turning.reset();
            continue;
        }

        const double requested_magnitude = std::min(
            std::fabs(proposed_step),
            config.max_mino_step);
        if (!std::isfinite(requested_magnitude) ||
            requested_magnitude < config.min_mino_step) {
            return terminate(
                current_public,
                rejection_reason(last_rejection),
                last_rejection == RejectionKind::None
                    ? "proposed Mino step is below the minimum"
                    : rejection_message(
                          last_rejection,
                          "fell below the minimum"));
        }
        const double attempted_step =
            integration_direction * requested_magnitude;
        if (current_mino + attempted_step == current_mino) {
            return terminate(
                current_public,
                TerminationReason::StepUnderflow,
                "proposed step cannot advance Mino time");
        }

        ++attempted_steps;
        if (turning_phase_authoritative) {
            const auto phase_step =
                detail::attempt_kerr_turning_phase_step(
                    *metric_,
                    constants,
                    potentials,
                    turning_phase,
                    current_mino,
                    attempted_step,
                    config);
            if (phase_step.status !=
                    numerics::Dopri5StepResult<7>::
                        Status::Completed ||
                !phase_step.accepted) {
                ++diagnostics.rejected_steps;
                ++consecutive_rejections;
                last_rejection =
                    phase_step.status ==
                            numerics::Dopri5StepResult<7>::
                                Status::Completed
                        ? RejectionKind::ErrorEstimate
                        : RejectionKind::NonFiniteStage;
                proposed_step =
                    phase_step.status ==
                            numerics::Dopri5StepResult<7>::
                                Status::Completed
                        ? phase_step.next_step
                        : attempted_step *
                              config.dopri5.min_factor;
                if (consecutive_rejections >=
                    config.max_rejections_per_step) {
                    return terminate(
                        current_public,
                        rejection_reason(last_rejection),
                        rejection_message(
                            last_rejection,
                            "reached the rejection limit"));
                }
                continue;
            }
            if (!phase_step.dense_output.has_value()) {
                return terminate(
                    current_public,
                    TerminationReason::NonFiniteState,
                    "accepted turning phase step has no dense output");
            }

            const auto selected_event =
                detail::select_first_kerr_phase_step_event(
                    *metric_,
                    constants,
                    current,
                    *phase_step.dense_output,
                    active_events);
            if (selected_event.status ==
                detail::KerrSeparatedEventStatus::Failed) {
                return terminate(
                    current_public,
                    TerminationReason::EventRootFailure,
                    selected_event.message);
            }

            const double accepted_mino =
                selected_event.hit.has_value()
                    ? selected_event.hit->mino_parameter
                    : current_mino + attempted_step;
            const auto accepted_phase =
                phase_step.dense_output->evaluate(
                    accepted_mino);
            const detail::KerrSeparatedState accepted =
                detail::project_kerr_turning_phase(
                    accepted_phase, current);
            PhaseSpaceState accepted_public;
            std::optional<EventHit> public_hit;
            if (selected_event.hit.has_value()) {
                accepted_public =
                    selected_event.hit->public_hit.state;
                public_hit =
                    selected_event.hit->public_hit;
            } else {
                try {
                    accepted_public =
                        detail::reconstruct_kerr_phase_space(
                            *metric_, constants, accepted);
                } catch (const std::exception& error) {
                    return terminate(
                        current_public,
                        TerminationReason::NonFiniteState,
                        std::string(
                            "accepted turning phase state failed: ") +
                            error.what());
                }
            }

            const double accepted_step_magnitude =
                std::fabs(accepted_mino - current_mino);
            if (public_hit.has_value() &&
                accepted_step_magnitude == 0.0) {
                return terminate(
                    accepted_public,
                    selected_event.reason,
                    "geodesic event is at the current state",
                    public_hit);
            }
            if (!finite_phase_space(accepted_public) ||
                !metric_->valid_point(accepted_public.x)) {
                return terminate(
                    current_public,
                    TerminationReason::InvalidMetricPoint,
                    "accepted turning phase state is outside Kerr BL");
            }
            double additional_minimum_radius =
                std::numeric_limits<double>::quiet_NaN();
            for (const detail::TurningCoordinate coordinate :
                 {detail::TurningCoordinate::Radial,
                  detail::TurningCoordinate::Polar}) {
                const bool locked =
                    coordinate ==
                            detail::TurningCoordinate::Radial
                        ? current.radial_direction ==
                              detail::SeparatedDirection::Locked
                        : current.polar_direction ==
                              detail::SeparatedDirection::Locked;
                if (locked) {
                    continue;
                }
                const auto crossing =
                    detail::locate_kerr_turning_phase_crossing(
                        coordinate,
                        potentials,
                        *phase_step.dense_output,
                        accepted_mino,
                        config.potential_tolerance,
                        config
                            .critical_derivative_tolerance);
                if (!crossing.has_value()) {
                    continue;
                }
                if (crossing->status ==
                    detail::TurningStatus::Failed) {
                    return terminate(
                        current_public,
                        TerminationReason::EventRootFailure,
                        crossing->message);
                }
                if (crossing->status ==
                    detail::TurningStatus::NearCritical) {
                    const auto critical_phase =
                        phase_step.dense_output->evaluate(
                            crossing->mino_parameter);
                    const auto critical_separated =
                        detail::project_kerr_turning_phase(
                            critical_phase, current);
                    PhaseSpaceState critical_public;
                    try {
                        critical_public =
                            detail::reconstruct_kerr_phase_space(
                                *metric_,
                                constants,
                                critical_separated);
                    } catch (const std::exception&) {
                        critical_public = current_public;
                    }
                    return terminate(
                        critical_public,
                        TerminationReason::NearCriticalOrbit,
                        crossing->message);
                }
                if (coordinate ==
                    detail::TurningCoordinate::Radial) {
                    ++diagnostics.radial_turns;
                    additional_minimum_radius =
                        crossing->root_radius_M;
                } else {
                    ++diagnostics.polar_turns;
                }
            }
            try {
                diagnostic_tracker.accept(
                    accepted,
                    accepted_public,
                    accepted_step_magnitude,
                    diagnostics,
                    additional_minimum_radius);
            } catch (const std::exception& error) {
                return terminate(
                    current_public,
                    TerminationReason::NonFiniteState,
                    std::string(
                        "turning phase diagnostics failed: ") +
                        error.what());
            }

            current = accepted;
            turning_phase = accepted_phase;
            current_public = accepted_public;
            current_mino = accepted_mino;
            if (public_hit.has_value()) {
                return terminate(
                    current_public,
                    selected_event.reason,
                    selected_event.message,
                    public_hit);
            }

            consecutive_rejections = 0;
            last_rejection = RejectionKind::None;
            proposed_step = phase_step.next_step;
            continue;
        }

        detail::KerrSeparatedStepAttempt trial =
            detail::attempt_kerr_separated_step(
                *metric_,
                constants,
                current,
                current_mino,
                attempted_step,
                config);
        const auto& step = trial.step;
        if (step.status !=
                numerics::Dopri5StepResult<5>::Status::Completed ||
            !step.accepted) {
            ++diagnostics.rejected_steps;
            ++consecutive_rejections;
            if (trial.radial_forbidden ||
                trial.polar_forbidden) {
                last_rejection =
                    RejectionKind::ForbiddenPotential;
                const detail::TurningCoordinate coordinate =
                    trial.radial_forbidden
                        ? detail::TurningCoordinate::Radial
                        : detail::TurningCoordinate::Polar;
                const double allowed_coordinate =
                    current.values[
                        coordinate ==
                                detail::TurningCoordinate::Radial
                            ? detail::kRadius
                            : detail::kMu];
                const double forbidden_coordinate =
                    coordinate ==
                            detail::TurningCoordinate::Radial
                        ? trial
                              .first_radial_forbidden_coordinate
                        : trial
                              .first_polar_forbidden_coordinate;
                const double fixed_other_coordinate =
                    current.values[
                        coordinate ==
                                detail::TurningCoordinate::Radial
                            ? detail::kMu
                            : detail::kRadius];
                const detail::TurningRoot root =
                    detail::locate_kerr_turning_root(
                        coordinate,
                        allowed_coordinate,
                        forbidden_coordinate,
                        potentials,
                        fixed_other_coordinate,
                        config.root_tolerance,
                        config.potential_tolerance,
                        config
                            .critical_derivative_tolerance);
                if (root.status ==
                    detail::TurningStatus::Failed) {
                    return terminate(
                        current_public,
                        TerminationReason::EventRootFailure,
                        root.message);
                }
                if (root.status ==
                    detail::TurningStatus::NearCritical) {
                    return terminate(
                        current_public,
                        TerminationReason::NearCriticalOrbit,
                        root.message);
                }
                pending_turning =
                    PendingTurning{coordinate, root};
                proposed_step =
                    attempted_step *
                    config.dopri5.min_factor;
            } else if (trial.invalid_metric_point) {
                last_rejection =
                    RejectionKind::InvalidMetricPoint;
                proposed_step =
                    attempted_step *
                    config.dopri5.min_factor;
            } else if (
                trial.non_finite_stage ||
                step.status !=
                    numerics::Dopri5StepResult<5>::
                        Status::Completed) {
                last_rejection =
                    RejectionKind::NonFiniteStage;
                proposed_step =
                    attempted_step *
                    config.dopri5.min_factor;
            } else {
                last_rejection =
                    RejectionKind::ErrorEstimate;
                proposed_step = step.next_step;
            }

            if (consecutive_rejections >=
                config.max_rejections_per_step) {
                return terminate(
                    current_public,
                    rejection_reason(last_rejection),
                    rejection_message(
                        last_rejection,
                        "reached the rejection limit"));
            }
            continue;
        }
        if (!step.dense_output.has_value()) {
            return terminate(
                current_public,
                TerminationReason::NonFiniteState,
                "accepted separated step has no dense output");
        }

        const detail::KerrSeparatedEventSelection selected_event =
            detail::select_first_kerr_step_event(
                *metric_,
                constants,
                current,
                *step.dense_output,
                active_events);
        if (selected_event.status ==
            detail::KerrSeparatedEventStatus::Failed) {
            return terminate(
                current_public,
                TerminationReason::EventRootFailure,
                selected_event.message);
        }

        detail::KerrSeparatedState accepted;
        PhaseSpaceState accepted_public;
        double accepted_mino;
        std::optional<EventHit> public_hit;
        if (selected_event.hit.has_value()) {
            accepted =
                selected_event.hit->separated_state;
            accepted_public =
                selected_event.hit->public_hit.state;
            accepted_mino =
                selected_event.hit->mino_parameter;
            public_hit =
                selected_event.hit->public_hit;
        } else {
            accepted = detail::KerrSeparatedState{
                step.state,
                current.radial_direction,
                current.polar_direction,
            };
            accepted_mino =
                current_mino + attempted_step;
            try {
                accepted_public =
                    detail::reconstruct_kerr_phase_space(
                        *metric_, constants, accepted);
            } catch (const std::exception& error) {
                return terminate(
                    current_public,
                    TerminationReason::NonFiniteState,
                    std::string(
                        "accepted separated state failed: ") +
                        error.what());
            }
        }

        const double accepted_step_magnitude =
            std::fabs(accepted_mino - current_mino);
        if (public_hit.has_value() &&
            accepted_step_magnitude == 0.0) {
            return terminate(
                accepted_public,
                selected_event.reason,
                "geodesic event is at the current state",
                public_hit);
        }
        if (!finite_phase_space(accepted_public) ||
            !metric_->valid_point(accepted_public.x)) {
            return terminate(
                current_public,
                TerminationReason::InvalidMetricPoint,
                "accepted separated state is outside Kerr BL");
        }

        try {
            turning_phase =
                detail::advance_kerr_turning_phase(
                    *metric_,
                    constants,
                    potentials,
                    turning_phase,
                    accepted_mino - current_mino,
                    config);
        } catch (const std::exception& error) {
            return terminate(
                current_public,
                TerminationReason::NonFiniteState,
                std::string(
                    "turning phase advance failed: ") +
                    error.what());
        }

        try {
            diagnostic_tracker.accept(
                accepted,
                accepted_public,
                accepted_step_magnitude,
                diagnostics);
        } catch (const std::exception& error) {
            return terminate(
                current_public,
                TerminationReason::NonFiniteState,
                std::string(
                    "accepted diagnostics failed: ") +
                    error.what());
        }

        current = accepted;
        current_public = accepted_public;
        current_mino = accepted_mino;

        if (public_hit.has_value()) {
            return terminate(
                current_public,
                selected_event.reason,
                selected_event.message,
                public_hit);
        }

        consecutive_rejections = 0;
        last_rejection = RejectionKind::None;
        proposed_step = step.next_step;
    }
}

} // namespace solar::relativity
