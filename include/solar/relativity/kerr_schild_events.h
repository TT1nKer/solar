#pragma once

#include "solar/relativity/geodesic_types.h"
#include "solar/relativity/kerr_schild_metric.h"

namespace solar::relativity {

double kerr_schild_interior_cutoff_radius(
    const KerrSchildCartesianMetric& metric,
    double configured_radius_M);

GeodesicEvent make_kerr_schild_horizon_event(
    const KerrSchildCartesianMetric& metric,
    double root_tolerance);

GeodesicEvent make_kerr_schild_interior_cutoff_event(
    const KerrSchildCartesianMetric& metric,
    double configured_radius_M,
    double root_tolerance);

} // namespace solar::relativity
