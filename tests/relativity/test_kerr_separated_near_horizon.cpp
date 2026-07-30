#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_separated.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/observer.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace solar::relativity;

namespace {

constexpr double half_pi =
    1.570796326794896619231321691639751442;
constexpr double pi =
    3.141592653589793238462643383279502884;

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << '\n';
    condition ? ++passed : ++failed;
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric metric(1.0, 0.0);
    const ObserverResult observer = make_zamo_observer(
        metric,
        Contravariant4{Vec4{{0.0, 30.0, half_pi, 0.0}}});
    check("Schwarzschild ZAMO is available", bool(observer));
    if (!observer) {
        return 1;
    }

    constexpr double pixel_x = 23.5;
    constexpr double pixel_y = 22.5;
    constexpr double resolution = 64.0;
    const double normalized_x =
        2.0 * pixel_x / resolution - 1.0;
    const double normalized_y =
        1.0 - 2.0 * pixel_y / resolution;
    const double tangent_half_fov =
        std::tan(20.0 * pi / 180.0);
    const InitialStateResult initialized =
        initialize_local_photon(
            metric,
            *observer.frame,
            Vec3{{
                1.0,
                normalized_y * tangent_half_fov,
                -normalized_x * tangent_half_fov,
            }});
    check("camera photon initializes", bool(initialized));
    if (!initialized) {
        return 1;
    }

    KerrSeparatedConfig config =
        KerrSeparatedConfig::cpu_reference(
            GeodesicKind::Null,
            1.0,
            -0.02,
            0.25,
            200.0);
    config.dopri5.relative_tolerance = 2.0e-13;
    config.potential_tolerance = 1.0e-8;
    config.root_tolerance = 1.0e-10;

    constexpr double capture_radius = 2.001;
    const GeodesicEvent capture{
        "consumer BL capture cutoff",
        [](const PhaseSpaceState& state) {
            return state.x.v[1] - capture_radius;
        },
        EventDirection::Decreasing,
        TerminationReason::InteriorCutoff,
        1.0e-10,
    };
    const GeodesicEvent escape{
        "consumer escape sphere",
        [](const PhaseSpaceState& state) {
            return state.x.v[1] - 60.0;
        },
        EventDirection::Increasing,
        TerminationReason::Escaped,
        1.0e-10,
    };
    const KerrSeparatedIntegrator integrator(metric);
    PhaseSpaceState current = *initialized.state;
    EventDirection equatorial_direction =
        EventDirection::Increasing;
    TerminationReason final_reason =
        TerminationReason::NonFiniteState;
    double max_constraint = 0.0;
    double max_carter = 0.0;
    double max_radial_residual = 0.0;
    std::size_t completed_segments = 0;
    for (std::size_t attempt = 0; attempt < 8; ++attempt) {
        const GeodesicEvent equatorial_return{
            "consumer equatorial return",
            [](const PhaseSpaceState& state) {
                return state.x.v[2] - half_pi;
            },
            equatorial_direction,
            TerminationReason::DiskSurfaceHit,
            1.0e-12,
        };
        const KerrSeparatedIntegrationResult result =
            integrator.integrate(
                current,
                config,
                {capture, escape, equatorial_return});
        ++completed_segments;
        max_constraint = std::max(
            max_constraint,
            result.diagnostics.max_constraint_error);
        max_carter = std::max(
            max_carter,
            result.diagnostics.max_carter_rel_error);
        max_radial_residual = std::max(
            max_radial_residual,
            result.diagnostics.max_radial_residual);
        current = result.final_state;
        final_reason = result.diagnostics.reason;
        if (final_reason != TerminationReason::DiskSurfaceHit) {
            break;
        }
        check(
            "equatorial event carries its exact state",
            result.event.has_value());
        if (!result.event) {
            return 1;
        }
        current = result.event->state;
        equatorial_direction =
            equatorial_direction == EventDirection::Increasing
                ? EventDirection::Decreasing
                : EventDirection::Increasing;
    }

    check(
        "consumer ray reaches the capture cutoff",
        final_reason ==
            TerminationReason::InteriorCutoff);
    check(
        "capture event returns its localized radius",
        std::fabs(current.x.v[1] - capture_radius) < 1.0e-9);
    check(
        "phase-event Hamiltonian diagnostic keeps headroom",
        max_constraint < 1.0e-11);
    check(
        "phase-event radial potential diagnostic keeps headroom",
        max_radial_residual < 1.0e-11);
    check(
        "phase-event Carter diagnostic stays below gate",
        max_carter < 1.0e-9);
    check(
        "observer-to-past affine parameter is negative",
        current.affine < 0.0);

    std::cout << std::setprecision(17)
              << "  segments=" << completed_segments
              << " final_radius=" << current.x.v[1]
              << " final_theta=" << current.x.v[2]
              << " final_pr=" << current.p.v[1]
              << " max_constraint=" << max_constraint
              << " max_carter=" << max_carter
              << " max_radial_residual="
              << max_radial_residual
              << " final_constraint="
              << hamiltonian_constraint_error(
                     metric, current, GeodesicKind::Null)
              << "\n\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
