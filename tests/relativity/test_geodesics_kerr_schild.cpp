#include "solar/constants.h"
#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_schild_events.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/observer.h"
#include "solar/relativity/spacetime_algebra.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(
    const char* name,
    bool condition,
    double diagnostic = 0.0) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << " (value=" << diagnostic << ")\n";
    condition ? ++passed : ++failed;
}

GeodesicIntegrationConfig plunge_config(
    double max_affine) {
    GeodesicIntegrationConfig config =
        GeodesicIntegrationConfig::cpu_reference(
            GeodesicKind::TimelikeUnitMass,
            1.0,
            2.0e-3,
            2.0e-2,
            max_affine);
    config.monitor_energy = true;
    config.monitor_lz = true;
    config.dopri5.relative_tolerance = 1.0e-12;
    for (double& tolerance :
         config.dopri5.absolute_tolerance) {
        tolerance *= 0.1;
    }
    config.stationary_energy_evaluator =
        kerr_schild_stationary_energy;
    config.axial_angular_momentum_evaluator =
        kerr_schild_axial_angular_momentum;
    return config;
}

bool finite_state(const PhaseSpaceState& state) {
    return std::isfinite(state.affine) &&
           state.x.v.all_finite() &&
           state.p.v.all_finite();
}

double maximum_state_difference(
    const PhaseSpaceState& left,
    const PhaseSpaceState& right) {
    double maximum = std::fabs(left.affine - right.affine);
    for (std::size_t component = 0; component < 4; ++component) {
        maximum = std::max(
            {maximum,
             std::fabs(
                 left.x.v[component] -
                 right.x.v[component]),
             std::fabs(
                 left.p.v[component] -
                 right.p.v[component])});
    }
    return maximum;
}

Contravariant4 normalized_timelike_tangent(
    const KerrSchildCartesianMetric& metric,
    const PhaseSpaceState& state) {
    Contravariant4 tangent = raise_index(
        metric.contravariant(state.x),
        state.p);
    const double norm = metric_inner_product(
        metric.covariant(state.x),
        tangent,
        tangent);
    tangent.v = tangent.v / std::sqrt(-norm);
    return tangent;
}

std::array<Contravariant4, 3> cartesian_spatial_seeds() {
    std::array<Contravariant4, 3> seeds{};
    seeds[0].v[1] = 1.0;
    seeds[1].v[2] = 1.0;
    seeds[2].v[3] = 1.0;
    return seeds;
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric bl_metric(1.0, 0.5);
    const KerrSchildCartesianMetric ks_metric(1.0, 0.5);
    const KerrChartTransform transform(1.0, 0.5);
    const Contravariant4 matching_position{
        Vec4{{
            0.0,
            6.0,
            0.5 * solar::constants::PI,
            0.0,
        }}};
    const ObserverResult zamo =
        make_zamo_observer(
            bl_metric, matching_position);
    check("safe matching ZAMO exists", bool(zamo));
    if (!zamo) {
        std::cout << "\n=== Results: " << passed
                  << " passed, " << ++failed
                  << " failed ===\n";
        return 1;
    }

    const InitialStateResult initialized =
        initialize_local_timelike(
            bl_metric,
            *zamo.frame,
            Vec3{{-0.75, 0.0, 0.0}});
    check("unit-mass BL plunge initializes",
          bool(initialized));
    if (!initialized) {
        std::cout << "\n=== Results: " << passed
                  << " passed, " << ++failed
                  << " failed ===\n";
        return 1;
    }
    const PhaseSpaceState initial_ks =
        transform.state_to_kerr_schild(
            *initialized.state);
    check("transformed plunge is finite",
          finite_state(initial_ks));
    check("transformed plunge constraint",
          hamiltonian_constraint_error(
              ks_metric,
              initial_ks,
              GeodesicKind::TimelikeUnitMass) < 1.0e-12,
          hamiltonian_constraint_error(
              ks_metric,
              initial_ks,
              GeodesicKind::TimelikeUnitMass));

    const double root_tolerance = 1.0e-10;
    const GeodesicIntegrationResult horizon =
        GeodesicIntegrator(ks_metric).integrate(
            initial_ks,
            plunge_config(30.0),
            {make_kerr_schild_horizon_event(
                ks_metric, root_tolerance)});
    check("plunge reaches outer horizon",
          horizon.diagnostics.reason ==
              TerminationReason::HorizonCrossing,
          static_cast<double>(horizon.diagnostics.reason));
    check("horizon event is retained",
          horizon.event.has_value());
    check("horizon state is finite",
          finite_state(horizon.final_state));
    const double horizon_radius =
        ks_metric.radial_coordinate(
            horizon.final_state.x);
    check("horizon radius localized",
          std::fabs(
              horizon_radius -
              ks_metric.outer_horizon_radius()) <
              5.0e-9,
          std::fabs(
              horizon_radius -
              ks_metric.outer_horizon_radius()));
    check("horizon Hamiltonian constraint gate",
          horizon.diagnostics.max_constraint_error <
              1.0e-10,
          horizon.diagnostics.max_constraint_error);
    check("horizon stationary energy drift gate",
          horizon.diagnostics.max_energy_rel_error <
              1.0e-12,
          horizon.diagnostics.max_energy_rel_error);
    check("horizon axial angular momentum drift gate",
          horizon.diagnostics.max_lz_rel_error <
              1.0e-12,
          horizon.diagnostics.max_lz_rel_error);
    check("horizon metric remains finite",
          all_finite(
              ks_metric.covariant(
                  horizon.final_state.x)) &&
              all_finite(
                  ks_metric.contravariant(
                      horizon.final_state.x)));

    const Contravariant4 horizon_velocity =
        normalized_timelike_tangent(
            ks_metric,
            horizon.final_state);
    const ObserverResult falling_observer =
        make_arbitrary_observer(
            ks_metric,
            horizon.final_state.x,
            horizon_velocity,
            cartesian_spatial_seeds());
    check("freely falling horizon observer exists",
          bool(falling_observer));
    const double tetrad_error =
        falling_observer
            ? tetrad_orthonormality_error(
                  ks_metric,
                  *falling_observer.frame)
            : std::numeric_limits<double>::infinity();
    check("horizon observer tetrad gate",
          tetrad_error < 1.0e-12,
          tetrad_error);

    const PhaseSpaceState restart_state =
        horizon.final_state;
    check("event state is exact restart state",
          horizon.event.has_value() &&
              maximum_state_difference(
                  restart_state,
                  horizon.event->state) == 0.0);
    const GeodesicIntegrationResult interior =
        GeodesicIntegrator(ks_metric).integrate(
            restart_state,
            plunge_config(30.0),
            {make_kerr_schild_interior_cutoff_event(
                ks_metric,
                0.0,
                root_tolerance)});
    check("continued plunge reaches interior cutoff",
          interior.diagnostics.reason ==
              TerminationReason::InteriorCutoff,
          static_cast<double>(interior.diagnostics.reason));
    check("interior event is retained",
          interior.event.has_value());
    check("interior state is finite",
          finite_state(interior.final_state));
    const double interior_radius =
        ks_metric.radial_coordinate(
            interior.final_state.x);
    check("interior radius localized at 0.05M",
          std::fabs(interior_radius - 0.05) < 5.0e-9,
          std::fabs(interior_radius - 0.05));
    check("proper/affine time increases through horizon",
          interior.final_state.affine >
              restart_state.affine,
          interior.final_state.affine -
              restart_state.affine);
    check("interior Hamiltonian constraint gate",
          interior.diagnostics.max_constraint_error <
              1.0e-10,
          interior.diagnostics.max_constraint_error);
    check("interior stationary energy drift gate",
          interior.diagnostics.max_energy_rel_error <
              1.0e-12,
          interior.diagnostics.max_energy_rel_error);
    check("interior axial angular momentum drift gate",
          interior.diagnostics.max_lz_rel_error <
              1.0e-12,
          interior.diagnostics.max_lz_rel_error);

    std::cout << "  horizon_affine="
              << horizon.final_state.affine
              << " horizon_radius=" << horizon_radius
              << " horizon_constraint="
              << horizon.diagnostics.max_constraint_error
              << " tetrad_error=" << tetrad_error
              << " interior_affine="
              << interior.final_state.affine
              << " interior_radius=" << interior_radius
              << " interior_constraint="
              << interior.diagnostics.max_constraint_error
              << '\n';
    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
