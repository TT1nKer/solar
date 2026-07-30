#include "solar/constants.h"
#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/observer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

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

double wrapped_difference(
    double left,
    double right) {
    return std::remainder(
        left - right,
        2.0 * solar::constants::PI);
}

double normalized_difference(
    double left,
    double right) {
    return std::fabs(left - right) /
           std::max(
               {1.0, std::fabs(left), std::fabs(right)});
}

double common_position_error(
    const Contravariant4& left,
    const Contravariant4& right,
    double mass) {
    return std::max(
        {std::fabs(left.v[0] - right.v[0]) / mass,
         std::fabs(left.v[1] - right.v[1]) / mass,
         std::fabs(left.v[2] - right.v[2]),
         std::fabs(wrapped_difference(
             left.v[3], right.v[3]))});
}

double common_momentum_error(
    const Covariant4& left,
    const Covariant4& right) {
    double maximum = 0.0;
    for (std::size_t component = 0; component < 4; ++component) {
        maximum = std::max(
            maximum,
            normalized_difference(
                left.v[component],
                right.v[component]));
    }
    return maximum;
}

double percentile_95(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t rank = static_cast<std::size_t>(
        std::ceil(0.95 * static_cast<double>(values.size())));
    return values[rank - 1];
}

struct CrosscheckFixture {
    double mass;
    double spin_chi;
    GeodesicKind kind;
    Contravariant4 initial_position;
    Vec3 local_motion;
    double target_radius;
    double max_affine;
};

struct CrosscheckResult {
    bool initialized = false;
    TerminationReason bl_reason =
        TerminationReason::NonFiniteState;
    TerminationReason ks_reason =
        TerminationReason::NonFiniteState;
    double position_error =
        std::numeric_limits<double>::infinity();
    double momentum_error =
        std::numeric_limits<double>::infinity();
    double max_constraint_error =
        std::numeric_limits<double>::infinity();
    double max_energy_error =
        std::numeric_limits<double>::infinity();
    double max_lz_error =
        std::numeric_limits<double>::infinity();
};

std::optional<PhaseSpaceState> initialize_fixture(
    const KerrBoyerLindquistMetric& metric,
    const CrosscheckFixture& fixture) {
    const ObserverResult observer =
        make_zamo_observer(
            metric, fixture.initial_position);
    if (!observer) {
        return std::nullopt;
    }
    const InitialStateResult initialized =
        fixture.kind == GeodesicKind::Null
            ? initialize_local_photon(
                  metric,
                  *observer.frame,
                  fixture.local_motion)
            : initialize_local_timelike(
                  metric,
                  *observer.frame,
                  fixture.local_motion);
    if (!initialized) {
        return std::nullopt;
    }
    return initialized.state;
}

GeodesicEvent bl_radius_event(double target_radius) {
    return GeodesicEvent{
        "BL common physical radius",
        [target_radius](const PhaseSpaceState& state) {
            return state.x.v[1] - target_radius;
        },
        EventDirection::Decreasing,
        TerminationReason::UserEvent,
        1.0e-10,
    };
}

GeodesicEvent ks_radius_event(
    const KerrSchildCartesianMetric& metric,
    double target_radius) {
    return GeodesicEvent{
        "KS common physical radius",
        [owned_metric = metric, target_radius](
            const PhaseSpaceState& state) {
            return owned_metric.radial_coordinate(state.x) -
                   target_radius;
        },
        EventDirection::Decreasing,
        TerminationReason::UserEvent,
        1.0e-10,
    };
}

CrosscheckResult run_crosscheck(
    const CrosscheckFixture& fixture) {
    const KerrBoyerLindquistMetric bl_metric(
        fixture.mass, fixture.spin_chi);
    const KerrSchildCartesianMetric ks_metric(
        fixture.mass, fixture.spin_chi);
    const KerrChartTransform transform(
        fixture.mass, fixture.spin_chi);
    const std::optional<PhaseSpaceState> bl_initial =
        initialize_fixture(bl_metric, fixture);
    if (!bl_initial.has_value()) {
        return {};
    }
    const PhaseSpaceState ks_initial =
        transform.state_to_kerr_schild(*bl_initial);

    GeodesicIntegrationConfig bl_config =
        GeodesicIntegrationConfig::cpu_reference(
            fixture.kind,
            fixture.mass,
            0.01 * fixture.mass,
            0.1 * fixture.mass,
            fixture.max_affine);
    bl_config.monitor_energy = true;
    bl_config.monitor_lz = true;
    GeodesicIntegrationConfig ks_config = bl_config;
    ks_config.stationary_energy_evaluator =
        kerr_schild_stationary_energy;
    ks_config.axial_angular_momentum_evaluator =
        kerr_schild_axial_angular_momentum;

    const GeodesicIntegrationResult bl =
        GeodesicIntegrator(bl_metric).integrate(
            *bl_initial,
            bl_config,
            {bl_radius_event(fixture.target_radius)});
    const GeodesicIntegrationResult ks =
        GeodesicIntegrator(ks_metric).integrate(
            ks_initial,
            ks_config,
            {ks_radius_event(
                ks_metric, fixture.target_radius)});

    CrosscheckResult result;
    result.initialized = true;
    result.bl_reason = bl.diagnostics.reason;
    result.ks_reason = ks.diagnostics.reason;
    result.max_constraint_error = std::max(
        bl.diagnostics.max_constraint_error,
        ks.diagnostics.max_constraint_error);
    result.max_energy_error = std::max(
        bl.diagnostics.max_energy_rel_error,
        ks.diagnostics.max_energy_rel_error);
    result.max_lz_error = std::max(
        bl.diagnostics.max_lz_rel_error,
        ks.diagnostics.max_lz_rel_error);
    if (bl.event.has_value() &&
        ks.event.has_value() &&
        bl.diagnostics.reason == TerminationReason::UserEvent &&
        ks.diagnostics.reason == TerminationReason::UserEvent) {
        const PhaseSpaceState ks_in_bl =
            transform.state_to_boyer_lindquist(
                ks.final_state);
        result.position_error = common_position_error(
            ks_in_bl.x,
            bl.final_state.x,
            fixture.mass);
        result.momentum_error = common_momentum_error(
            ks_in_bl.p,
            bl.final_state.p);
    }
    return result;
}

} // namespace

int main() {
    const std::vector<CrosscheckFixture> ordinary{
        {1.0, 0.5, GeodesicKind::Null,
         Contravariant4{
             Vec4{{0.0, 10.0, 1.1, 0.2}}},
         Vec3{{-1.0, 0.10, 0.15}}, 8.0, 20.0},
        {1.0, -0.5, GeodesicKind::Null,
         Contravariant4{
             Vec4{{1.0, 9.0, 2.0, -0.4}}},
         Vec3{{-1.0, -0.12, -0.20}}, 7.0, 20.0},
        {1.5, 0.0, GeodesicKind::Null,
         Contravariant4{
             Vec4{{-2.0, 12.0, 1.3, 1.0}}},
         Vec3{{-1.0, 0.05, 0.10}}, 9.5, 30.0},
        {1.0, 0.4, GeodesicKind::TimelikeUnitMass,
         Contravariant4{
             Vec4{{0.0, 10.0, 1.0, 0.5}}},
         Vec3{{-0.35, 0.03, 0.08}}, 9.0, 30.0},
        {1.0, -0.6, GeodesicKind::TimelikeUnitMass,
         Contravariant4{
             Vec4{{0.5, 11.0, 1.8, -0.8}}},
         Vec3{{-0.40, -0.04, -0.10}}, 9.5, 30.0},
        {2.0, 0.7, GeodesicKind::TimelikeUnitMass,
         Contravariant4{
             Vec4{{-1.0, 24.0, 1.3, 2.0}}},
         Vec3{{-0.45, 0.05, 0.15}}, 21.0, 60.0},
    };

    std::vector<double> position_errors;
    std::vector<double> momentum_errors;
    double maximum_constraint_error = 0.0;
    double maximum_energy_error = 0.0;
    double maximum_lz_error = 0.0;
    for (const CrosscheckFixture& fixture : ordinary) {
        const CrosscheckResult result =
            run_crosscheck(fixture);
        check("ordinary fixture initializes",
              result.initialized);
        check("ordinary BL reaches common event",
              result.bl_reason == TerminationReason::UserEvent,
              static_cast<double>(result.bl_reason));
        check("ordinary KS reaches common event",
              result.ks_reason == TerminationReason::UserEvent,
              static_cast<double>(result.ks_reason));
        check("ordinary termination classifications agree",
              result.bl_reason == result.ks_reason);
        position_errors.push_back(result.position_error);
        momentum_errors.push_back(result.momentum_error);
        maximum_constraint_error = std::max(
            maximum_constraint_error,
            result.max_constraint_error);
        maximum_energy_error = std::max(
            maximum_energy_error,
            result.max_energy_error);
        maximum_lz_error = std::max(
            maximum_lz_error,
            result.max_lz_error);
    }

    const double position_p95 =
        percentile_95(position_errors);
    const double momentum_p95 =
        percentile_95(momentum_errors);
    check("ordinary BL/KS position P95 gate",
          position_p95 < 1.0e-8,
          position_p95);
    check("ordinary BL/KS momentum P95 gate",
          momentum_p95 < 1.0e-8,
          momentum_p95);
    check("ordinary Hamiltonian constraint gate",
          maximum_constraint_error < 1.0e-10,
          maximum_constraint_error);
    check("ordinary stationary energy drift gate",
          maximum_energy_error < 1.0e-12,
          maximum_energy_error);
    check("ordinary axial angular momentum drift gate",
          maximum_lz_error < 1.0e-12,
          maximum_lz_error);

    const KerrSchildCartesianMetric near_metric(1.0, 0.5);
    const CrosscheckFixture near_horizon{
        1.0,
        0.5,
        GeodesicKind::Null,
        Contravariant4{
            Vec4{{0.0, 3.0, 1.2, 0.3}}},
        Vec3{{-1.0, 0.02, 0.05}},
        near_metric.outer_horizon_radius() + 0.02,
        20.0,
    };
    const CrosscheckResult near =
        run_crosscheck(near_horizon);
    check("near-horizon fixture initializes",
          near.initialized);
    check("near-horizon BL reaches common event",
          near.bl_reason == TerminationReason::UserEvent,
          static_cast<double>(near.bl_reason));
    check("near-horizon KS reaches common event",
          near.ks_reason == TerminationReason::UserEvent,
          static_cast<double>(near.ks_reason));
    check("near-horizon termination classifications agree",
          near.bl_reason == near.ks_reason);
    check("near-horizon position gate",
          near.position_error < 1.0e-6,
          near.position_error);
    check("near-horizon momentum gate",
          near.momentum_error < 1.0e-6,
          near.momentum_error);
    check("near-horizon constraint gate",
          near.max_constraint_error < 1.0e-10,
          near.max_constraint_error);

    std::cout << "  ordinary_position_p95=" << position_p95
              << " ordinary_momentum_p95=" << momentum_p95
              << " max_constraint=" << maximum_constraint_error
              << " max_energy=" << maximum_energy_error
              << " max_lz=" << maximum_lz_error
              << " near_position=" << near.position_error
              << " near_momentum=" << near.momentum_error
              << '\n';
    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
