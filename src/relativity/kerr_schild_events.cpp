#include "solar/relativity/kerr_schild_events.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

void require_root_tolerance(double root_tolerance) {
    if (!std::isfinite(root_tolerance) ||
        root_tolerance <= 0.0) {
        throw std::invalid_argument(
            "Kerr-Schild event root tolerance must be "
            "finite and positive");
    }
}

GeodesicEvent make_decreasing_radius_event(
    const char* name,
    const KerrSchildCartesianMetric& metric,
    double target_radius,
    TerminationReason reason,
    double root_tolerance) {
    require_root_tolerance(root_tolerance);
    return GeodesicEvent{
        name,
        [owned_metric = metric, target_radius](
            const PhaseSpaceState& state) {
            return owned_metric.radial_coordinate(state.x) -
                   target_radius;
        },
        EventDirection::Decreasing,
        reason,
        root_tolerance,
    };
}

} // namespace

double kerr_schild_interior_cutoff_radius(
    const KerrSchildCartesianMetric& metric,
    double configured_radius_M) {
    if (!std::isfinite(configured_radius_M) ||
        configured_radius_M < 0.0) {
        throw std::invalid_argument(
            "Kerr-Schild interior cutoff radius must be "
            "finite and non-negative");
    }
    return std::max(
        0.05 * metric.mass(),
        configured_radius_M);
}

GeodesicEvent make_kerr_schild_horizon_event(
    const KerrSchildCartesianMetric& metric,
    double root_tolerance) {
    return make_decreasing_radius_event(
        "Kerr-Schild outer horizon",
        metric,
        metric.outer_horizon_radius(),
        TerminationReason::HorizonCrossing,
        root_tolerance);
}

GeodesicEvent make_kerr_schild_interior_cutoff_event(
    const KerrSchildCartesianMetric& metric,
    double configured_radius_M,
    double root_tolerance) {
    return make_decreasing_radius_event(
        "Kerr-Schild interior cutoff",
        metric,
        kerr_schild_interior_cutoff_radius(
            metric, configured_radius_M),
        TerminationReason::InteriorCutoff,
        root_tolerance);
}

} // namespace solar::relativity
