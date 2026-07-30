#pragma once

#include "solar/relativity/metric.h"

namespace solar::relativity::detail {

struct ThinDiskSurfaceGeometry {
    double radius = 0.0;
    Contravariant4 normal;
};

ThinDiskSurfaceGeometry evaluate_thin_disk_surface_geometry(
    const Metric& metric,
    const Contravariant4& position);

bool valid_thin_disk_surface_geometry(
    const Metric& metric,
    const Contravariant4& position,
    const Contravariant4& normal,
    const Contravariant4& emitter_four_velocity);

} // namespace solar::relativity::detail
