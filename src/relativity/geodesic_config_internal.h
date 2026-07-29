#pragma once

#include "solar/relativity/geodesic_integrator.h"

namespace solar::relativity::detail {

void validate_geodesic_config(
    const GeodesicIntegrationConfig& config);

} // namespace solar::relativity::detail
