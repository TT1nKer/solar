#pragma once

#include "solar/relativity/kerr_bl_metric.h"

#include <cstddef>
#include <vector>

namespace solar::relativity {

struct ShadowCriticalPoint {
    double alpha;
    double beta;
    double photon_radius;
};

std::vector<ShadowCriticalPoint> bardeen_shadow_curve(
    const KerrBoyerLindquistMetric& metric,
    double inclination,
    std::size_t samples_per_branch);

} // namespace solar::relativity
