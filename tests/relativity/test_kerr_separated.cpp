#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/kerr_separated.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace solar::relativity;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  PASS: " << name << "\n";
        ++passed;
    } else {
        std::cerr << "  FAIL: " << name << "\n";
        ++failed;
    }
}

void check_near(
    const std::string& name,
    double actual,
    double expected,
    double tolerance) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

PhaseSpaceState radial_photon(
    const KerrBoyerLindquistMetric& metric,
    double radius,
    double radial_sign) {
    PhaseSpaceState state{
        0.0,
        Contravariant4{
            Vec4{{0.0, radius, kPi / 2.0, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0, 0.0}}},
    };
    const Mat4 inverse = metric.contravariant(state.x);
    const double radial_squared =
        -(inverse[0][0] *
          state.p.v[0] * state.p.v[0]) /
        inverse[1][1];
    state.p.v[1] =
        std::copysign(std::sqrt(radial_squared), radial_sign);
    return state;
}

PhaseSpaceState constrained_photon(
    const KerrBoyerLindquistMetric& metric,
    double radius,
    double theta,
    double lz,
    double p_theta,
    double radial_sign) {
    PhaseSpaceState state{
        0.0,
        Contravariant4{
            Vec4{{0.0, radius, theta, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, p_theta, lz}}},
    };
    const Mat4 inverse = metric.contravariant(state.x);
    double non_radial_twice_hamiltonian = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        if (row == 1) {
            continue;
        }
        for (std::size_t column = 0; column < 4; ++column) {
            if (column == 1) {
                continue;
            }
            non_radial_twice_hamiltonian +=
                inverse[row][column] *
                state.p.v[row] *
                state.p.v[column];
        }
    }
    const double radial_squared =
        -non_radial_twice_hamiltonian /
        inverse[1][1];
    if (!(radial_squared >= 0.0) ||
        !std::isfinite(radial_squared)) {
        throw std::domain_error(
            "test photon has no real radial momentum");
    }
    state.p.v[1] =
        std::copysign(std::sqrt(radial_squared), radial_sign);
    return state;
}

KerrSeparatedConfig reference_config(
    double initial_step = 1.0e-4,
    double max_step = 1.0e-3,
    double max_affine = 20.0) {
    return KerrSeparatedConfig::cpu_reference(
        GeodesicKind::Null,
        1.0,
        initial_step,
        max_step,
        max_affine);
}

GeodesicEvent radius_event(
    const char* name,
    double radius,
    TerminationReason reason) {
    return GeodesicEvent{
        name,
        [radius](const PhaseSpaceState& state) {
            return state.x.v[1] - radius;
        },
        EventDirection::Increasing,
        reason,
        1.0e-10,
    };
}

template <typename Function>
void check_invalid_argument(
    const std::string& name,
    Function&& function) {
    try {
        function();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    }
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric schwarzschild(1.0, 0.0);
    const KerrSeparatedIntegrator integrator(schwarzschild);
    const PhaseSpaceState outgoing =
        radial_photon(schwarzschild, 10.0, 1.0);
    check(
        "radial photon fixture satisfies Hamiltonian constraint",
        hamiltonian_constraint_error(
            schwarzschild,
            outgoing,
            GeodesicKind::Null) <
            1.0e-14);

    const auto escaped = integrator.integrate(
        outgoing,
        reference_config(),
        {radius_event(
            "escape sphere", 12.0, TerminationReason::Escaped)});
    check(
        "escape event reason",
        escaped.diagnostics.reason ==
            TerminationReason::Escaped);
    check("escape event exists", escaped.event.has_value());
    check_near(
        "escape radius",
        escaped.final_state.x.v[1],
        12.0,
        1.0e-9);
    check(
        "ordinary ray accepts steps",
        escaped.diagnostics.accepted_steps > 0);
    check(
        "ordinary ray has no radial turn",
        escaped.diagnostics.radial_turns == 0);
    check(
        "equatorial plane remains locked",
        std::fabs(
            escaped.final_state.x.v[2] - kPi / 2.0) <
            1.0e-13);

    const GeodesicEvent initial_event{
        "initial radius",
        [](const PhaseSpaceState& state) {
            return state.x.v[1] - 10.0;
        },
        EventDirection::Any,
        TerminationReason::UserEvent,
        1.0e-10,
    };
    const auto initial_hit = integrator.integrate(
        outgoing, reference_config(), {initial_event});
    check(
        "initial any event terminates",
        initial_hit.diagnostics.reason ==
            TerminationReason::UserEvent);
    check(
        "initial any event accepts no step",
        initial_hit.diagnostics.accepted_steps == 0);
    check(
        "initial any event rejects no step",
        initial_hit.diagnostics.rejected_steps == 0);

    const auto first_event = integrator.integrate(
        outgoing,
        reference_config(5.0e-3, 5.0e-3),
        {
            radius_event(
                "near sphere",
                10.2,
                TerminationReason::MaterialSurfaceHit),
            radius_event(
                "far sphere",
                10.4,
                TerminationReason::Escaped),
        });
    check(
        "first event in integration direction wins",
        first_event.diagnostics.reason ==
            TerminationReason::MaterialSurfaceHit);
    check(
        "first event index retained",
        first_event.event.has_value() &&
            first_event.event->event_index == 0);
    check_near(
        "first event radius",
        first_event.final_state.x.v[1],
        10.2,
        1.0e-9);

    const auto affine_limit = integrator.integrate(
        outgoing, reference_config(1.0e-4, 1.0e-3, 0.1));
    check(
        "affine limit reason",
        affine_limit.diagnostics.reason ==
            TerminationReason::MaxAffine);
    check_near(
        "affine limit displacement",
        std::fabs(
            affine_limit.final_state.affine -
            outgoing.affine),
        0.1,
        1.0e-10);

    auto coordinate_config =
        reference_config(1.0e-4, 1.0e-3, 20.0);
    coordinate_config.max_coordinate_time = 0.05;
    const auto coordinate_limit = integrator.integrate(
        outgoing, coordinate_config);
    check(
        "coordinate-time limit reason",
        coordinate_limit.diagnostics.reason ==
            TerminationReason::MaxCoordinateTime);
    check_near(
        "coordinate-time displacement",
        std::fabs(
            coordinate_limit.final_state.x.v[0] -
            outgoing.x.v[0]),
        0.05,
        1.0e-9);

    auto step_config =
        reference_config(1.0e-5, 1.0e-5, 20.0);
    step_config.max_total_steps = 1;
    const auto step_limit =
        integrator.integrate(outgoing, step_config);
    check(
        "total-step limit reason",
        step_limit.diagnostics.reason ==
            TerminationReason::MaxSteps);
    check(
        "total-step limit accepts one step",
        step_limit.diagnostics.accepted_steps == 1);

    const auto backward = integrator.integrate(
        outgoing,
        reference_config(-1.0e-4, 1.0e-3, 0.1));
    check(
        "negative Mino integration reaches affine limit",
        backward.diagnostics.reason ==
            TerminationReason::MaxAffine);
    check(
        "negative Mino integration decreases affine",
        backward.final_state.affine < outgoing.affine);
    check(
        "negative Mino integration follows past radial branch",
        backward.final_state.x.v[1] <
            outgoing.x.v[1]);
    check(
        "backward-traced momentum remains future directed",
        backward.final_state.p.v[0] < 0.0);

    const PhaseSpaceState scattering =
        constrained_photon(
            schwarzschild,
            10.0,
            kPi / 2.0,
            6.0,
            0.0,
            -1.0);
    const auto scattered = integrator.integrate(
        scattering,
        reference_config(1.0e-4, 1.0e-3, 100.0),
        {radius_event(
            "return sphere",
            10.0,
            TerminationReason::Escaped)});
    check(
        "scattering ray returns through outer sphere",
        scattered.diagnostics.reason ==
            TerminationReason::Escaped);
    check(
        "scattering ray counts one radial turn",
        scattered.diagnostics.radial_turns == 1);
    check(
        "scattering minimum radius is finite and interior",
        std::isfinite(
            scattered.diagnostics.min_radius_M) &&
            scattered.diagnostics.min_radius_M < 5.0 &&
            scattered.diagnostics.min_radius_M > 4.0);

    const auto scattered_fine = integrator.integrate(
        scattering,
        reference_config(5.0e-5, 5.0e-4, 100.0),
        {radius_event(
            "fine return sphere",
            10.0,
            TerminationReason::Escaped)});
    check(
        "fine scattering ray also returns",
        scattered_fine.diagnostics.reason ==
            TerminationReason::Escaped);
    check(
        "turning return converges under halved step",
        std::fabs(
            scattered.final_state.x.v[0] -
            scattered_fine.final_state.x.v[0]) <
            1.0e-7 &&
            std::fabs(
                scattered.final_state.x.v[3] -
                scattered_fine.final_state.x.v[3]) <
            1.0e-7);

    const double starting_mu = 0.5;
    const double starting_theta = std::acos(starting_mu);
    const PhaseSpaceState polar_turn =
        constrained_photon(
            schwarzschild,
            10.0,
            starting_theta,
            2.0,
            -0.1,
            1.0);
    const GeodesicEvent polar_return{
        "polar return",
        [starting_theta](const PhaseSpaceState& state) {
            return state.x.v[2] - starting_theta;
        },
        EventDirection::Increasing,
        TerminationReason::UserEvent,
        1.0e-10,
    };
    const auto polar = integrator.integrate(
        polar_turn,
        reference_config(1.0e-5, 2.0e-4, 100.0),
        {polar_return});
    check(
        "polar return event reached",
        polar.diagnostics.reason ==
            TerminationReason::UserEvent);
    check(
        "polar ray counts one polar turn",
        polar.diagnostics.polar_turns == 1);
    check_near(
        "polar ray returns to starting theta",
        polar.final_state.x.v[2],
        starting_theta,
        1.0e-9);

    const PhaseSpaceState critical_root{
        0.0,
        Contravariant4{
            Vec4{{0.0, 3.0, kPi / 2.0, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0, std::sqrt(27.0)}}},
    };
    const auto critical = integrator.integrate(
        critical_root,
        reference_config(1.0e-5, 1.0e-4, 10.0));
    check(
        "spherical photon is explicitly critical",
        critical.diagnostics.reason ==
            TerminationReason::NearCriticalOrbit);
    check(
        "critical photon is not a simple radial turn",
        critical.diagnostics.radial_turns == 0);
    std::cout
        << std::setprecision(17)
        << "  turning_probe radial_reason="
        << static_cast<int>(scattered.diagnostics.reason)
        << " radial_r=" << scattered.final_state.x.v[1]
        << " radial_min=" << scattered.diagnostics.min_radius_M
        << " radial_rejected="
        << scattered.diagnostics.rejected_steps
        << " radial_accepted="
        << scattered.diagnostics.accepted_steps
        << " radial_message="
        << scattered.diagnostics.message
        << " convergence_dt="
        << std::fabs(
               scattered.final_state.x.v[0] -
               scattered_fine.final_state.x.v[0])
        << " convergence_dphi="
        << std::fabs(
               scattered.final_state.x.v[3] -
               scattered_fine.final_state.x.v[3])
        << " polar_reason="
        << static_cast<int>(polar.diagnostics.reason)
        << " polar_theta=" << polar.final_state.x.v[2]
        << " polar_rejected="
        << polar.diagnostics.rejected_steps
        << " polar_accepted="
        << polar.diagnostics.accepted_steps
        << " polar_message="
        << polar.diagnostics.message << "\n";

    const GeodesicEvent malformed_event{
        "missing callback",
        {},
        EventDirection::Increasing,
        TerminationReason::UserEvent,
        1.0e-10,
    };
    const auto malformed = integrator.integrate(
        outgoing, reference_config(), {malformed_event});
    check(
        "malformed event is explicit",
        malformed.diagnostics.reason ==
            TerminationReason::EventRootFailure);
    check(
        "malformed event performs no integration",
        malformed.diagnostics.accepted_steps == 0 &&
            malformed.diagnostics.rejected_steps == 0);

    GeodesicEvent unknown_direction = initial_event;
    unknown_direction.direction =
        static_cast<EventDirection>(91);
    const auto invalid_direction = integrator.integrate(
        outgoing, reference_config(), {unknown_direction});
    check(
        "unknown event direction is explicit",
        invalid_direction.diagnostics.reason ==
            TerminationReason::EventRootFailure);
    check(
        "unknown event direction performs no integration",
        invalid_direction.diagnostics.accepted_steps == 0 &&
            invalid_direction.diagnostics.rejected_steps == 0);

    check_invalid_argument(
        "unknown geodesic kind rejected",
        [&] {
            auto invalid = reference_config();
            invalid.kind = static_cast<GeodesicKind>(91);
            static_cast<void>(
                integrator.integrate(outgoing, invalid));
        });
    check_invalid_argument(
        "zero initial Mino step rejected",
        [&] {
            static_cast<void>(
                KerrSeparatedConfig::cpu_reference(
                    GeodesicKind::Null,
                    1.0,
                    0.0,
                    1.0e-3,
                    1.0));
        });
    check_invalid_argument(
        "negative maximum Mino step rejected",
        [&] {
            static_cast<void>(
                KerrSeparatedConfig::cpu_reference(
                    GeodesicKind::Null,
                    1.0,
                    1.0e-4,
                    -1.0e-3,
                    1.0));
        });
    check_invalid_argument(
        "non-finite potential tolerance rejected",
        [&] {
            auto invalid = reference_config();
            invalid.potential_tolerance =
                std::numeric_limits<double>::quiet_NaN();
            static_cast<void>(
                integrator.integrate(outgoing, invalid));
        });

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
