#include "solar/relativity/kerr_shadow.h"
#include "solar/relativity/kerr_separated.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/observer.h"
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
              << "}\n";
    return 0;
}
