#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(const char* name, double actual, double expected,
                double tolerance) {
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) <= tolerance);
}

template <typename Action>
void check_invalid_argument(const char* name, Action action) {
    try {
        action();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    }
}

PhaseSpaceState minkowski_photon() {
    return PhaseSpaceState{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
}

PhaseSpaceState minkowski_massive() {
    return PhaseSpaceState{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.25, 0.75, 0.0, 0.0}}},
    };
}

GeodesicIntegrationConfig null_config(
    double initial_step,
    double max_step,
    double max_affine) {
    return GeodesicIntegrationConfig::cpu_reference(
        GeodesicKind::Null,
        1.0,
        initial_step,
        max_step,
        max_affine);
}

} // namespace

int main() {
    const MinkowskiMetric metric;
    const GeodesicIntegrator integrator(metric);
    const PhaseSpaceState photon = minkowski_photon();

    const auto null_line = integrator.integrate(
        photon, null_config(0.1, 0.5, 10.0));
    check("null line reaches affine limit",
          null_line.diagnostics.reason ==
              TerminationReason::MaxAffine);
    check_near("null affine limit",
               null_line.final_state.affine, 10.0, 0.0);
    check_near("null coordinate time",
               null_line.final_state.x.v[0], 10.0, 2.0e-14);
    check_near("null spatial line",
               null_line.final_state.x.v[1], 10.0, 2.0e-14);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near("null momentum remains constant",
                   null_line.final_state.p.v[component],
                   photon.p.v[component], 0.0);
    }
    check("null constraint remains exact",
          null_line.diagnostics.max_constraint_error < 1.0e-14);
    check("null integration accepted steps",
          null_line.diagnostics.accepted_steps > 0);
    check("null integration has no rejected steps",
          null_line.diagnostics.rejected_steps == 0);

    auto timelike_config =
        GeodesicIntegrationConfig::cpu_reference(
            GeodesicKind::TimelikeUnitMass,
            1.0,
            0.2,
            1.0,
            10.0);
    timelike_config.max_proper_time = 4.0;
    const auto timelike_line = integrator.integrate(
        minkowski_massive(), timelike_config);
    check("timelike line reaches proper-time limit",
          timelike_line.diagnostics.reason ==
              TerminationReason::MaxProperTime);
    check_near("timelike affine proper time",
               timelike_line.final_state.affine, 4.0, 2.0e-12);
    check_near("timelike coordinate time",
               timelike_line.final_state.x.v[0], 5.0, 2.0e-12);
    check_near("timelike spatial line",
               timelike_line.final_state.x.v[1], 3.0, 2.0e-12);
    check_near("timelike final Hamiltonian",
               hamiltonian(metric, timelike_line.final_state),
               -0.5, 2.0e-14);

    const auto forward = integrator.integrate(
        photon, null_config(0.2, 0.7, 3.0));
    const auto backward = integrator.integrate(
        forward.final_state, null_config(-0.2, 0.7, 3.0));
    check("backward line reaches affine limit",
          backward.diagnostics.reason ==
              TerminationReason::MaxAffine);
    check_near("forward-backward affine",
               backward.final_state.affine, photon.affine, 2.0e-12);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near("forward-backward coordinate",
                   backward.final_state.x.v[component],
                   photon.x.v[component], 2.0e-12);
        check_near("forward-backward momentum",
                   backward.final_state.p.v[component],
                   photon.p.v[component], 2.0e-12);
    }

    const GeodesicEvent user_event{
        "x=2.25",
        [](const PhaseSpaceState& state) {
            return state.x.v[1] - 2.25;
        },
        EventDirection::Increasing,
        TerminationReason::UserEvent,
        1.0e-12,
    };
    const auto event_result = integrator.integrate(
        photon, null_config(3.0, 3.0, 10.0), {user_event});
    check("user event terminates with declared reason",
          event_result.diagnostics.reason ==
              TerminationReason::UserEvent);
    check("user event payload returned",
          event_result.event.has_value());
    check_near("user event affine",
               event_result.final_state.affine, 2.25, 1.0e-12);
    check_near("user event coordinate",
               event_result.final_state.x.v[1], 2.25, 1.0e-12);

    const GeodesicEvent initial_event{
        "initial x",
        [](const PhaseSpaceState& state) {
            return state.x.v[1];
        },
        EventDirection::Any,
        TerminationReason::UserEvent,
        1.0e-12,
    };
    const auto initial_event_result = integrator.integrate(
        photon, null_config(1.0, 1.0, 10.0), {initial_event});
    check("initial event terminates",
          initial_event_result.diagnostics.reason ==
              TerminationReason::UserEvent);
    check("initial event accepts no step",
          initial_event_result.diagnostics.accepted_steps == 0);
    check("initial event keeps min step unavailable",
          std::isnan(initial_event_result.diagnostics.min_step));

    const std::vector<GeodesicEvent> competing_events{
        GeodesicEvent{
            "later",
            [](const PhaseSpaceState& state) {
                return state.x.v[1] - 1.5;
            },
            EventDirection::Increasing,
            TerminationReason::DiskSurfaceHit,
            1.0e-12,
        },
        GeodesicEvent{
            "earlier",
            [](const PhaseSpaceState& state) {
                return state.x.v[1] - 1.0;
            },
            EventDirection::Increasing,
            TerminationReason::MaterialSurfaceHit,
            1.0e-12,
        },
    };
    const auto first_event = integrator.integrate(
        photon,
        null_config(2.0, 2.0, 10.0),
        competing_events);
    check("first event in step wins",
          first_event.diagnostics.reason ==
              TerminationReason::MaterialSurfaceHit);
    check("first event index retained",
          first_event.event.has_value() &&
              first_event.event->event_index == 1);
    check_near("first event affine",
               first_event.final_state.affine, 1.0, 1.0e-12);

    auto coordinate_limit = null_config(2.0, 2.0, 10.0);
    coordinate_limit.max_coordinate_time = 0.75;
    const auto coordinate_limited =
        integrator.integrate(photon, coordinate_limit);
    check("coordinate-time limit reason",
          coordinate_limited.diagnostics.reason ==
              TerminationReason::MaxCoordinateTime);
    check_near("coordinate-time limit state",
               coordinate_limited.final_state.x.v[0],
               0.75, 1.0e-12);

    auto step_limit = null_config(0.1, 0.1, 10.0);
    step_limit.max_total_steps = 1;
    const auto step_limited =
        integrator.integrate(photon, step_limit);
    check("total-step limit reason",
          step_limited.diagnostics.reason ==
              TerminationReason::MaxSteps);
    check("total-step limit accepted one step",
          step_limited.diagnostics.accepted_steps == 1);

    auto underflow_config = null_config(0.1, 1.0, 10.0);
    underflow_config.min_step = 0.2;
    const auto underflow =
        integrator.integrate(photon, underflow_config);
    check("underflowing initial step is explicit",
          underflow.diagnostics.reason ==
              TerminationReason::StepUnderflow);
    check("underflow accepts no steps",
          underflow.diagnostics.accepted_steps == 0);

    PhaseSpaceState constraint_violation = photon;
    constraint_violation.p.v[1] = 0.0;
    const auto invalid_constraint = integrator.integrate(
        constraint_violation, null_config(0.1, 1.0, 1.0));
    check("initial constraint violation is explicit",
          invalid_constraint.diagnostics.reason ==
              TerminationReason::ConstraintViolation);
    check("constraint violation accepts no steps",
          invalid_constraint.diagnostics.accepted_steps == 0);

    PhaseSpaceState non_finite = photon;
    non_finite.x.v[2] =
        std::numeric_limits<double>::quiet_NaN();
    const auto invalid_state = integrator.integrate(
        non_finite, null_config(0.1, 1.0, 1.0));
    check("non-finite initial state is explicit",
          invalid_state.diagnostics.reason ==
              TerminationReason::NonFiniteState);

    const SchwarzschildBoyerLindquistMetric schwarzschild(1.0);
    const GeodesicIntegrator schwarzschild_integrator(schwarzschild);
    PhaseSpaceState invalid_point{
        0.0,
        Contravariant4{
            Vec4{{0.0, 2.0, 1.5707963267948966, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
    const auto invalid_domain = schwarzschild_integrator.integrate(
        invalid_point, null_config(0.1, 1.0, 1.0));
    check("invalid initial metric point is explicit",
          invalid_domain.diagnostics.reason ==
              TerminationReason::InvalidMetricPoint);

    auto zero_step_config = null_config(0.1, 1.0, 1.0);
    zero_step_config.initial_step = 0.0;
    check_invalid_argument("zero initial step rejected", [&] {
        (void)integrator.integrate(photon, zero_step_config);
    });

    auto null_proper_limit = null_config(0.1, 1.0, 1.0);
    null_proper_limit.max_proper_time = 0.5;
    check_invalid_argument("null proper-time limit rejected", [&] {
        (void)integrator.integrate(photon, null_proper_limit);
    });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
