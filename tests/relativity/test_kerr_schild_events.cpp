#include "solar/relativity/kerr_schild_events.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(
    const char* name,
    double actual,
    double expected,
    double tolerance) {
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) <= tolerance);
}

PhaseSpaceState state_at_radius(
    const KerrSchildCartesianMetric& metric,
    double radius) {
    return PhaseSpaceState{
        0.0,
        Contravariant4{
            Vec4{{0.0, radius, metric.spin_length(), 0.0}}},
        Covariant4{},
    };
}

template <typename Callable>
void check_invalid_argument(
    const char* name,
    Callable&& callable) {
    try {
        callable();
        check(name, false);
    } catch (const std::invalid_argument&) {
        check(name, true);
    }
}

GeodesicEvent event_from_temporary_metric() {
    return make_kerr_schild_horizon_event(
        KerrSchildCartesianMetric(1.0, 0.5),
        1.0e-11);
}

} // namespace

int main() {
    const KerrSchildCartesianMetric metric(2.0, 0.6);
    check_near(
        "zero configuration selects default cutoff",
        kerr_schild_interior_cutoff_radius(metric, 0.0),
        0.1,
        0.0);
    check_near(
        "configured radius below floor selects default",
        kerr_schild_interior_cutoff_radius(metric, 0.01),
        0.1,
        0.0);
    check_near(
        "configured radius above floor is retained",
        kerr_schild_interior_cutoff_radius(metric, 0.25),
        0.25,
        0.0);

    const GeodesicEvent horizon =
        make_kerr_schild_horizon_event(metric, 1.0e-11);
    check("horizon event is named", !horizon.name.empty());
    check("horizon event is decreasing",
          horizon.direction == EventDirection::Decreasing);
    check("horizon event reason",
          horizon.reason ==
              TerminationReason::HorizonCrossing);
    check_near(
        "horizon event root",
        horizon.function(
            state_at_radius(
                metric,
                metric.outer_horizon_radius())),
        0.0,
        2.0e-15);
    check("horizon event is positive outside",
          horizon.function(
              state_at_radius(
                  metric,
                  metric.outer_horizon_radius() + 0.2)) > 0.0);
    check("horizon event is negative inside",
          horizon.function(
              state_at_radius(
                  metric,
                  metric.outer_horizon_radius() - 0.2)) < 0.0);

    const GeodesicEvent interior =
        make_kerr_schild_interior_cutoff_event(
            metric, 0.0, 2.0e-12);
    check("interior event is named", !interior.name.empty());
    check("interior event is decreasing",
          interior.direction == EventDirection::Decreasing);
    check("interior event reason",
          interior.reason ==
              TerminationReason::InteriorCutoff);
    check_near(
        "interior event root",
        interior.function(state_at_radius(metric, 0.1)),
        0.0,
        2.0e-15);
    check_near(
        "event root tolerance is retained",
        interior.root_tolerance,
        2.0e-12,
        0.0);

    const GeodesicEvent owned =
        event_from_temporary_metric();
    const KerrSchildCartesianMetric matching_metric(1.0, 0.5);
    check_near(
        "event owns temporary metric value",
        owned.function(
            state_at_radius(
                matching_metric,
                matching_metric.outer_horizon_radius())),
        0.0,
        2.0e-15);

    check_invalid_argument(
        "negative configured cutoff rejected",
        [&] {
            (void)kerr_schild_interior_cutoff_radius(
                metric, -0.1);
        });
    check_invalid_argument(
        "non-finite configured cutoff rejected",
        [&] {
            (void)kerr_schild_interior_cutoff_radius(
                metric,
                std::numeric_limits<double>::infinity());
        });
    check_invalid_argument(
        "zero horizon root tolerance rejected",
        [&] {
            (void)make_kerr_schild_horizon_event(
                metric, 0.0);
        });
    check_invalid_argument(
        "non-finite interior root tolerance rejected",
        [&] {
            (void)make_kerr_schild_interior_cutoff_event(
                metric,
                0.0,
                std::numeric_limits<double>::quiet_NaN());
        });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
