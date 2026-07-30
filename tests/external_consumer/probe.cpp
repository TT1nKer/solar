#include "solar/relativity/emission_model.h"
#include "solar/relativity/fluid_model.h"
#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/kerr_shadow.h"
#include "solar/relativity/kerr_chart_transform.h"
#include "solar/relativity/kerr_schild_events.h"
#include "solar/relativity/kerr_schild_metric.h"
#include "solar/relativity/kerr_separated.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/observer.h"
#include "solar/relativity/radiative_transfer.h"
#include "solar/relativity/thin_disk.h"
#include "solar/version.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    constexpr double half_pi = 1.5707963267948966;
    constexpr std::size_t samples_per_branch = 65;
    constexpr double expected_left_edge = -4.096266658713869;
    constexpr double expected_right_edge = 6.138155724715452;
    constexpr double edge_tolerance = 1.0e-13;

    if (solar::version != "0.2.0-alpha.1" ||
        solar::physics_contract != "relativity-v3-phase2") {
        std::cerr << "unexpected Solar public contract\n";
        return 1;
    }

    const solar::relativity::KerrBoyerLindquistMetric metric(1.0, 0.5);
    const auto curve = solar::relativity::bardeen_shadow_curve(
        metric, half_pi, samples_per_branch);
    if (curve.size() != 2 * samples_per_branch - 2) {
        std::cerr << "unexpected Kerr shadow sample count\n";
        return 2;
    }

    const auto extrema = std::minmax_element(
        curve.begin(),
        curve.end(),
        [](const auto& left, const auto& right) {
            return left.alpha < right.alpha;
        });
    const double left_edge = extrema.first->alpha;
    const double right_edge = extrema.second->alpha;
    if (!std::isfinite(left_edge) ||
        !std::isfinite(right_edge) ||
        std::abs(left_edge - expected_left_edge) > edge_tolerance ||
        std::abs(right_edge - expected_right_edge) > edge_tolerance) {
        std::cerr << "unexpected Kerr shadow edge\n";
        return 3;
    }

    const solar::relativity::Contravariant4 observer_position{
        solar::relativity::Vec4{{
            0.0,
            20.0,
            half_pi,
            0.0,
        }}};
    const auto observer =
        solar::relativity::make_zamo_observer(
            metric, observer_position);
    if (!observer) {
        std::cerr << "installed ZAMO construction failed\n";
        return 4;
    }
    const auto photon =
        solar::relativity::initialize_local_photon(
            metric,
            *observer.frame,
            solar::relativity::Vec3{{-1.0, 0.0, 0.0}});
    if (!photon) {
        std::cerr << "installed photon initialization failed\n";
        return 5;
    }
    const auto separated_config =
        solar::relativity::KerrSeparatedConfig::cpu_reference(
            solar::relativity::GeodesicKind::Null,
            1.0,
            1.0e-5,
            1.0e-4,
            0.1);
    const auto separated =
        solar::relativity::KerrSeparatedIntegrator(metric)
            .integrate(*photon.state, separated_config);
    if (separated.diagnostics.reason !=
            solar::relativity::TerminationReason::MaxAffine ||
        separated.diagnostics.accepted_steps == 0 ||
        separated.diagnostics.max_constraint_error >= 1.0e-10) {
        std::cerr << "installed separated Kerr integration failed\n";
        return 6;
    }

    const solar::relativity::KerrSchildCartesianMetric
        kerr_schild_metric(1.0, 0.5);
    const solar::relativity::KerrChartTransform transform(
        1.0, 0.5);
    const solar::relativity::PhaseSpaceState kerr_schild_photon =
        transform.state_to_kerr_schild(*photon.state);
    const auto covariant =
        kerr_schild_metric.covariant(kerr_schild_photon.x);
    const auto contravariant =
        kerr_schild_metric.contravariant(
            kerr_schild_photon.x);
    const auto identity =
        solar::relativity::multiply(
            covariant, contravariant);
    double identity_error = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected =
                row == column ? 1.0 : 0.0;
            identity_error = std::max(
                identity_error,
                std::abs(identity[row][column] - expected));
        }
    }
    if (kerr_schild_metric.chart() !=
            solar::relativity::Chart::KerrSchildCartesian ||
        identity_error >= 5.0e-13) {
        std::cerr << "installed Kerr-Schild metric failed\n";
        return 7;
    }

    auto kerr_schild_config =
        solar::relativity::GeodesicIntegrationConfig::cpu_reference(
            solar::relativity::GeodesicKind::Null,
            1.0,
            1.0e-3,
            1.0e-2,
            5.0e-2);
    kerr_schild_config.monitor_energy = true;
    kerr_schild_config.monitor_lz = true;
    kerr_schild_config.stationary_energy_evaluator =
        solar::relativity::kerr_schild_stationary_energy;
    kerr_schild_config.axial_angular_momentum_evaluator =
        solar::relativity::kerr_schild_axial_angular_momentum;
    const auto kerr_schild =
        solar::relativity::GeodesicIntegrator(
            kerr_schild_metric)
            .integrate(
                kerr_schild_photon,
                kerr_schild_config);
    const auto horizon_event =
        solar::relativity::make_kerr_schild_horizon_event(
            kerr_schild_metric, 1.0e-10);
    if (kerr_schild.diagnostics.reason !=
            solar::relativity::TerminationReason::MaxAffine ||
        kerr_schild.diagnostics.accepted_steps == 0 ||
        kerr_schild.diagnostics.max_constraint_error >= 1.0e-10 ||
        kerr_schild.diagnostics.max_energy_rel_error >= 1.0e-12 ||
        kerr_schild.diagnostics.max_lz_rel_error >= 1.0e-12 ||
        horizon_event.direction !=
            solar::relativity::EventDirection::Decreasing ||
        horizon_event.function(kerr_schild_photon) <= 0.0) {
        std::cerr << "installed Kerr-Schild integration failed\n";
        return 8;
    }

    const solar::relativity::MinkowskiMetric minkowski;
    const solar::relativity::PhaseSpaceState flat_photon{
        0.0,
        solar::relativity::Contravariant4{
            solar::relativity::Vec4{{0.0, 0.0, 0.0, 0.0}}},
        solar::relativity::Covariant4{
            solar::relativity::Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
    const solar::relativity::Contravariant4 static_emitter{
        solar::relativity::Vec4{{1.0, 0.0, 0.0, 0.0}}};
    const auto flat_redshift =
        solar::relativity::evaluate_redshift(
            minkowski,
            flat_photon,
            static_emitter,
            1.0);
    const auto transfer =
        solar::relativity::advance_backward_transfer(
            {},
            solar::relativity::TransferCoefficients{2.0, 0.5},
            3.0);
    constexpr double expected_transfer_intensity =
        3.1074793594062806;
    constexpr double expected_transfer_transmission =
        0.22313016014842982;

    const solar::relativity::Contravariant4 disk_point{
        solar::relativity::Vec4{{
            0.0, 8.0, half_pi, 0.2}}};
    const solar::relativity::AnalyticCircularDiskConfig
        disk_config{
            1.0,
            0.5,
            solar::relativity::OrbitSense::Prograde,
            6.0,
            20.0,
            1.0,
            16.407349347422414,
            0.0,
            1.0e-8,
        };
    const solar::relativity::AnalyticCircularDiskFluid
        disk(disk_config);
    const auto disk_sample = disk.sample(metric, disk_point);
    const solar::relativity::AnalyticOpticallyThinTorus torus(
        solar::relativity::AnalyticOpticallyThinTorusConfig{
            1.0,
            0.5,
            solar::relativity::OrbitSense::Prograde,
            8.0,
            2.0,
            0.2,
            3.0,
            4.0,
            0.0,
            1.0e-4,
        });
    const auto torus_sample =
        torus.sample(metric, disk_point);
    const auto disk_emitter =
        solar::relativity::make_equatorial_circular_observer(
            metric,
            disk_point,
            solar::relativity::OrbitSense::Prograde);
    if (!disk_emitter) {
        std::cerr << "installed disk emitter construction failed\n";
        return 9;
    }
    const auto disk_photon =
        solar::relativity::initialize_local_photon(
            metric,
            *disk_emitter.frame,
            solar::relativity::Vec3{{0.0, -1.0, 0.0}});
    if (!disk_photon) {
        std::cerr << "installed disk photon construction failed\n";
        return 10;
    }
    solar::relativity::ThinDiskCrossingRecorder surface(
        solar::relativity::ThinDiskRecorderConfig{
            solar::relativity::DiskOpacityMode::Opaque, 8},
        disk,
        solar::relativity::ThinDiskSurfaceEmission(
            0.75, 10.0 / 4096.0, 0.7));
    const auto surface_result =
        surface.record(metric, *disk_photon.state, 1.0);

    const double transfer_error = std::max(
        std::abs(
            transfer.state.invariant_intensity -
            expected_transfer_intensity),
        std::abs(
            transfer.state.transmission -
            expected_transfer_transmission));
    if (!flat_redshift ||
        std::abs(flat_redshift.sample.redshift_g - 1.0) >
            5.0e-14 ||
        !transfer ||
        !std::isfinite(transfer.state.invariant_intensity) ||
        !std::isfinite(transfer.state.transmission) ||
        transfer_error >= 5.0e-14 ||
        !disk_sample.valid ||
        !torus_sample.valid ||
        !surface_result ||
        !surface_result.recorded ||
        surface.crossings().size() != 1 ||
        !std::isfinite(disk_sample.temperature) ||
        !std::isfinite(torus_sample.density) ||
        !std::isfinite(
            surface.observed().specific_intensity)) {
        std::cerr << "installed Phase 5 transfer API failed\n";
        return 11;
    }

    std::cout << std::setprecision(17)
              << "{\"solar_version\":\"" << solar::version
              << "\",\"physics_contract\":\"" << solar::physics_contract
              << "\",\"samples\":" << curve.size()
              << ",\"left_edge\":" << left_edge
              << ",\"right_edge\":" << right_edge
              << ",\"separated_steps\":"
              << separated.diagnostics.accepted_steps
              << ",\"separated_constraint\":"
              << separated.diagnostics.max_constraint_error
              << ",\"ks_steps\":"
              << kerr_schild.diagnostics.accepted_steps
              << ",\"ks_constraint\":"
              << kerr_schild.diagnostics.max_constraint_error
              << ",\"ks_inverse_error\":"
              << identity_error
              << ",\"transfer_intensity\":"
              << transfer.state.invariant_intensity
              << ",\"transfer_transmission\":"
              << transfer.state.transmission
              << ",\"disk_temperature\":"
              << disk_sample.temperature
              << ",\"torus_density\":"
              << torus_sample.density
              << ",\"surface_specific_intensity\":"
              << surface.observed().specific_intensity
              << ",\"surface_crossings\":"
              << surface.crossings().size()
              << "}\n";
    return 0;
}
