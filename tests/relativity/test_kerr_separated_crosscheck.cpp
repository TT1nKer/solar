#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/kerr_constants.h"
#include "solar/relativity/kerr_separated.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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

PhaseSpaceState constrained_state(
    const KerrBoyerLindquistMetric& metric,
    GeodesicKind kind,
    double radius,
    double theta,
    double energy,
    double lz,
    double p_theta,
    double radial_sign) {
    PhaseSpaceState state{
        0.0,
        Contravariant4{
            Vec4{{0.0, radius, theta, 0.0}}},
        Covariant4{
            Vec4{{-energy, 0.0, p_theta, lz}}},
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
        (2.0 * hamiltonian_target(kind) -
         non_radial_twice_hamiltonian) /
        inverse[1][1];
    if (!(radial_squared >= 0.0) ||
        !std::isfinite(radial_squared)) {
        throw std::domain_error(
            "cross-check fixture has no real radial momentum");
    }
    state.p.v[1] =
        std::copysign(std::sqrt(radial_squared), radial_sign);
    return state;
}

GeodesicIntegrationConfig hamiltonian_config(
    const KerrBoyerLindquistMetric& metric,
    GeodesicKind kind,
    double max_affine) {
    auto config =
        GeodesicIntegrationConfig::cpu_reference(
            kind, metric.mass(), 1.0e-3, 2.0e-2, max_affine);
    config.dopri5.relative_tolerance = 2.0e-13;
    config.constraint_tolerance = 1.0e-10;
    config.carter_evaluator =
        [&metric, kind](const PhaseSpaceState& state) {
            return evaluate_kerr_constants(
                       metric, state, kind)
                .Q;
        };
    return config;
}

KerrSeparatedConfig separated_config(
    const KerrBoyerLindquistMetric& metric,
    GeodesicKind kind,
    double max_mino_step,
    double max_affine) {
    auto config = KerrSeparatedConfig::cpu_reference(
        kind,
        metric.mass(),
        0.1 * max_mino_step,
        max_mino_step,
        max_affine);
    config.dopri5.relative_tolerance = 2.0e-13;
    return config;
}

double coordinate_error_M(
    const PhaseSpaceState& separated,
    const PhaseSpaceState& hamiltonian,
    double mass_M) {
    const double dt =
        (separated.x.v[0] - hamiltonian.x.v[0]) /
        mass_M;
    const double dr =
        (separated.x.v[1] - hamiltonian.x.v[1]) /
        mass_M;
    const double dtheta =
        separated.x.v[2] - hamiltonian.x.v[2];
    const double dphi = std::remainder(
        separated.x.v[3] - hamiltonian.x.v[3],
        2.0 * kPi);
    return std::sqrt(
        dt * dt + dr * dr +
        dtheta * dtheta + dphi * dphi);
}

double carter_relative_error(
    const KerrBoyerLindquistMetric& metric,
    GeodesicKind kind,
    const PhaseSpaceState& initial,
    const PhaseSpaceState& final) {
    const double initial_q =
        evaluate_kerr_constants(metric, initial, kind).Q;
    const double final_q =
        evaluate_kerr_constants(metric, final, kind).Q;
    return std::fabs(final_q - initial_q) /
           std::max(
               std::fabs(initial_q),
               metric.mass() * metric.mass() * 1.0e-14);
}

struct FixtureResult {
    std::string name;
    double coarse_error;
    double fine_error;
    double affine_error;
    double separated_constraint;
    double hamiltonian_constraint;
    double separated_carter;
    double hamiltonian_carter;
    double delta_t;
    double delta_r;
    double delta_theta;
    double delta_phi;
    std::size_t radial_turns;
    std::size_t polar_turns;
    TerminationReason separated_reason;
    TerminationReason hamiltonian_reason;
};

FixtureResult run_fixture(
    std::string name,
    const KerrBoyerLindquistMetric& metric,
    GeodesicKind kind,
    const PhaseSpaceState& initial,
    const GeodesicEvent& event,
    double coarse_mino_step,
    double max_affine) {
    const GeodesicIntegrator hamiltonian_integrator(metric);
    const auto hamiltonian = hamiltonian_integrator.integrate(
        initial,
        hamiltonian_config(metric, kind, max_affine),
        {event});

    const KerrSeparatedIntegrator separated_integrator(metric);
    const auto coarse = separated_integrator.integrate(
        initial,
        separated_config(
            metric, kind, coarse_mino_step, max_affine),
        {event});
    const auto fine = separated_integrator.integrate(
        initial,
        separated_config(
            metric,
            kind,
            0.5 * coarse_mino_step,
            max_affine),
        {event});
    return FixtureResult{
        std::move(name),
        coordinate_error_M(
            coarse.final_state,
            hamiltonian.final_state,
            metric.mass()),
        coordinate_error_M(
            fine.final_state,
            hamiltonian.final_state,
            metric.mass()),
        std::fabs(
            coarse.final_state.affine -
            hamiltonian.final_state.affine) /
            metric.mass(),
        coarse.diagnostics.max_constraint_error,
        hamiltonian.diagnostics.max_constraint_error,
        coarse.diagnostics.max_carter_rel_error,
        carter_relative_error(
            metric,
            kind,
            initial,
            hamiltonian.final_state),
        coarse.final_state.x.v[0] -
            hamiltonian.final_state.x.v[0],
        coarse.final_state.x.v[1] -
            hamiltonian.final_state.x.v[1],
        coarse.final_state.x.v[2] -
            hamiltonian.final_state.x.v[2],
        std::remainder(
            coarse.final_state.x.v[3] -
                hamiltonian.final_state.x.v[3],
            2.0 * kPi),
        coarse.diagnostics.radial_turns,
        coarse.diagnostics.polar_turns,
        coarse.diagnostics.reason,
        hamiltonian.diagnostics.reason,
    };
}

GeodesicEvent increasing_radius_event(
    const char* name,
    double radius) {
    return GeodesicEvent{
        name,
        [radius](const PhaseSpaceState& state) {
            return state.x.v[1] - radius;
        },
        EventDirection::Increasing,
        TerminationReason::UserEvent,
        1.0e-11,
    };
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric prograde_metric(1.0, 0.5);
    const KerrBoyerLindquistMetric retrograde_metric(1.0, -0.5);
    const KerrBoyerLindquistMetric schwarzschild(1.0, 0.0);

    std::vector<FixtureResult> fixtures;
    fixtures.push_back(run_fixture(
        "ordinary-null-positive-spin",
        prograde_metric,
        GeodesicKind::Null,
        constrained_state(
            prograde_metric,
            GeodesicKind::Null,
            12.0,
            1.1,
            1.0,
            2.0,
            0.7,
            1.0),
        increasing_radius_event("r=14", 14.0),
        2.0e-4,
        50.0));
    fixtures.push_back(run_fixture(
        "ordinary-null-negative-spin",
        retrograde_metric,
        GeodesicKind::Null,
        constrained_state(
            retrograde_metric,
            GeodesicKind::Null,
            12.0,
            1.1,
            1.0,
            2.0,
            0.7,
            1.0),
        increasing_radius_event("r=14", 14.0),
        2.0e-4,
        50.0));
    fixtures.push_back(run_fixture(
        "ordinary-timelike-positive-spin",
        prograde_metric,
        GeodesicKind::TimelikeUnitMass,
        constrained_state(
            prograde_metric,
            GeodesicKind::TimelikeUnitMass,
            12.0,
            1.2,
            1.05,
            2.0,
            0.5,
            1.0),
        increasing_radius_event("r=13", 13.0),
        2.0e-4,
        50.0));
    fixtures.push_back(run_fixture(
        "ordinary-null-schwarzschild",
        schwarzschild,
        GeodesicKind::Null,
        constrained_state(
            schwarzschild,
            GeodesicKind::Null,
            10.0,
            1.2,
            1.0,
            1.0,
            0.5,
            1.0),
        increasing_radius_event("r=12", 12.0),
        2.0e-4,
        50.0));

    const PhaseSpaceState scattering = constrained_state(
        schwarzschild,
        GeodesicKind::Null,
        10.0,
        kPi / 2.0,
        1.0,
        6.0,
        0.0,
        -1.0);
    fixtures.push_back(run_fixture(
        "radial-turn-schwarzschild",
        schwarzschild,
        GeodesicKind::Null,
        scattering,
        increasing_radius_event("return-r=10", 10.0),
        2.0e-4,
        100.0));

    const double polar_theta = std::acos(0.5);
    const PhaseSpaceState polar = constrained_state(
        schwarzschild,
        GeodesicKind::Null,
        10.0,
        polar_theta,
        1.0,
        2.0,
        -0.1,
        1.0);
    const GeodesicEvent polar_return{
        "polar return",
        [polar_theta](const PhaseSpaceState& state) {
            return state.x.v[2] - polar_theta;
        },
        EventDirection::Increasing,
        TerminationReason::UserEvent,
        1.0e-11,
    };
    fixtures.push_back(run_fixture(
        "polar-turn-schwarzschild",
        schwarzschild,
        GeodesicKind::Null,
        polar,
        polar_return,
        1.0e-4,
        100.0));

    std::vector<double> errors;
    double maximum_error = 0.0;
    for (const FixtureResult& fixture : fixtures) {
        check(
            fixture.name + " separated reaches common event",
            fixture.separated_reason ==
                TerminationReason::UserEvent);
        check(
            fixture.name + " Hamiltonian reaches common event",
            fixture.hamiltonian_reason ==
                TerminationReason::UserEvent);
        check(
            fixture.name + " common-event coordinate error",
            fixture.coarse_error < 1.0e-8);
        check(
            fixture.name + " affine reparameterization recovery",
            fixture.affine_error < 1.0e-8);
        check(
            fixture.name + " separated constraint gate",
            fixture.separated_constraint < 1.0e-10);
        check(
            fixture.name + " Hamiltonian constraint gate",
            fixture.hamiltonian_constraint < 1.0e-10);
        check(
            fixture.name + " separated Carter gate",
            fixture.separated_carter < 1.0e-10);
        check(
            fixture.name + " Hamiltonian Carter gate",
            fixture.hamiltonian_carter < 1.0e-10);
        check(
            fixture.name + " halved-step convergence",
            fixture.fine_error <=
                1.1 * fixture.coarse_error + 1.0e-11);
        errors.push_back(fixture.coarse_error);
        maximum_error = std::max(
            maximum_error, fixture.coarse_error);
        std::cout
            << "  fixture=" << fixture.name
            << " coarse_error=" << fixture.coarse_error
            << " fine_error=" << fixture.fine_error
            << " affine_error=" << fixture.affine_error
            << " separated_constraint="
            << fixture.separated_constraint
            << " hamiltonian_constraint="
            << fixture.hamiltonian_constraint
            << " separated_carter="
            << fixture.separated_carter
            << " hamiltonian_carter="
            << fixture.hamiltonian_carter
            << " dt=" << fixture.delta_t
            << " dr=" << fixture.delta_r
            << " dtheta=" << fixture.delta_theta
            << " dphi=" << fixture.delta_phi
            << " radial_turns=" << fixture.radial_turns
            << " polar_turns=" << fixture.polar_turns
            << "\n";
    }

    std::sort(errors.begin(), errors.end());
    const std::size_t p95_index = static_cast<std::size_t>(
        std::ceil(0.95 * static_cast<double>(errors.size()))) - 1;
    const double p95_error = errors[p95_index];
    check("worldline P95 gate", p95_error < 1.0e-8);
    check("worldline maximum gate", maximum_error < 1.0e-7);
    check(
        "radial fixture records one radial turn",
        fixtures[4].radial_turns == 1);
    check(
        "polar fixture records one polar turn",
        fixtures[5].polar_turns == 1);

    std::cout << "  worldline_p95=" << p95_error
              << " worldline_max=" << maximum_error << "\n";
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
