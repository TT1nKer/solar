#include "solar/relativity/event_root.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace solar::numerics;
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

Dopri5Config<8> exact_config() {
    Dopri5Config<8> config{};
    config.absolute_tolerance.fill(1.0e-12);
    config.relative_tolerance = 1.0e-11;
    config.safety = 0.9;
    config.min_factor = 0.2;
    config.max_factor = 5.0;
    config.error_norm = ErrorNorm::RootMeanSquare;
    return config;
}

GeodesicEvent radial_event(
    double root,
    EventDirection direction) {
    return GeodesicEvent{
        "radial-test",
        [root](const PhaseSpaceState& state) {
            return state.x.v[1] - root;
        },
        direction,
        TerminationReason::UserEvent,
        1.0e-12,
    };
}

} // namespace

int main() {
    StateN<8> initial{};
    const auto constant_rhs = [](double, const StateN<8>&) {
        StateN<8> derivative{};
        derivative[1] = 1.0;
        return derivative;
    };
    const auto step = dopri5_step(
        initial, 0.0, 1.0, constant_rhs, exact_config());
    check("constant dense step accepted", step.accepted);
    check("constant dense output available",
          step.dense_output.has_value());
    const auto dense = step.dense_output.value();

    const auto increasing = locate_event(
        4, radial_event(0.3, EventDirection::Increasing), dense);
    check("increasing event found",
          increasing.status == EventRootStatus::Found);
    check("increasing event carries hit", increasing.hit.has_value());
    check_near("increasing event affine",
               increasing.hit->affine, 0.3, 1.0e-12);
    check_near("increasing event state",
               increasing.hit->state.x.v[1], 0.3, 1.0e-12);
    check("event index retained",
          increasing.hit->event_index == 4);

    const auto wrong_direction = locate_event(
        0, radial_event(0.3, EventDirection::Decreasing), dense);
    check("direction mismatch is no root",
          wrong_direction.status == EventRootStatus::NoRoot);
    check("no-root result has no hit",
          !wrong_direction.hit.has_value());

    const auto unknown_direction = locate_event(
        0,
        radial_event(
            0.3, static_cast<EventDirection>(99)),
        dense);
    check("unknown event direction fails",
          unknown_direction.status == EventRootStatus::Failed);

    const auto any_direction = locate_event(
        1, radial_event(0.3, EventDirection::Any), dense);
    check("any-direction event found",
          any_direction.status == EventRootStatus::Found);
    check_near("any-direction affine",
               any_direction.hit->affine, 0.3, 1.0e-12);

    const auto start_root = locate_event(
        2, radial_event(0.0, EventDirection::Any), dense);
    check("start endpoint root found",
          start_root.status == EventRootStatus::Found);
    check_near("start endpoint affine",
               start_root.hit->affine, 0.0, 0.0);

    const auto end_root = locate_event(
        3, radial_event(1.0, EventDirection::Any), dense);
    check("end endpoint root found",
          end_root.status == EventRootStatus::Found);
    check_near("end endpoint affine",
               end_root.hit->affine, 1.0, 0.0);

    StateN<8> backward_initial{};
    backward_initial[1] = 1.0;
    const auto backward_rhs = [](double, const StateN<8>&) {
        StateN<8> derivative{};
        derivative[1] = 1.0;
        return derivative;
    };
    const auto backward_step = dopri5_step(
        backward_initial, 0.0, -1.0,
        backward_rhs, exact_config());
    const auto backward = locate_event(
        5, radial_event(0.3, EventDirection::Decreasing),
        backward_step.dense_output.value());
    check("negative-step decreasing event found",
          backward.status == EventRootStatus::Found);
    check_near("negative-step event affine",
               backward.hit->affine, -0.7, 1.0e-12);
    check_near("negative-step event state",
               backward.hit->state.x.v[1], 0.3, 1.0e-12);

    GeodesicEvent non_finite_internal = radial_event(
        0.3, EventDirection::Any);
    non_finite_internal.function =
        [](const PhaseSpaceState& state) {
            if (state.x.v[1] > 0.1 && state.x.v[1] < 0.9) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return state.x.v[1] - 0.3;
        };
    const auto failed_root = locate_event(
        6, non_finite_internal, dense);
    check("non-finite internal event value fails",
          failed_root.status == EventRootStatus::Failed);
    check("failed root has no hit",
          !failed_root.hit.has_value());

    GeodesicEvent unresolvable_root{
        "unresolvable tolerance",
        [](const PhaseSpaceState& state) {
            const double difference = state.x.v[1] - 0.3;
            if (std::fabs(difference) < 1.0e-14) {
                return std::copysign(
                    std::numeric_limits<double>::denorm_min(),
                    difference == 0.0 ? 1.0 : difference);
            }
            return difference;
        },
        EventDirection::Any,
        TerminationReason::UserEvent,
        std::numeric_limits<double>::denorm_min(),
    };
    const auto exhausted_root = locate_event(
        7, unresolvable_root, dense);
    check("unresolvable root tolerance fails explicitly",
          exhausted_root.status == EventRootStatus::Failed);
    check("exhausted root has no hit",
          !exhausted_root.hit.has_value());

    const IntegrationDiagnostics diagnostics{};
    check("diagnostic min step starts unavailable",
          std::isnan(diagnostics.min_step));
    check("diagnostic max step starts unavailable",
          std::isnan(diagnostics.max_step));
    check("energy diagnostic starts unavailable",
          std::isnan(diagnostics.max_energy_rel_error));
    check("Lz diagnostic starts unavailable",
          std::isnan(diagnostics.max_lz_rel_error));
    check("Carter diagnostic starts unavailable",
          std::isnan(diagnostics.max_carter_rel_error));
    check("Carter absolute diagnostic starts unavailable",
          std::isnan(diagnostics.max_carter_abs_error));
    check("v3 interior cutoff reason is representable",
          TerminationReason::InteriorCutoff !=
              TerminationReason::HorizonCrossing);

    std::cout << std::setprecision(17)
              << "  linear_event_affine_error="
              << std::fabs(increasing.hit->affine - 0.3)
              << " negative_step_event_affine_error="
              << std::fabs(backward.hit->affine + 0.7)
              << '\n';
    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
