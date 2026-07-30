#pragma once

#include "solar/relativity/fluid_model.h"

namespace solar::relativity::detail {

struct KerrFluidPoint {
    Contravariant4 boyer_lindquist_position;
    Contravariant4 caller_four_velocity;
    double radius;
    double theta;
    double equatorial_height;
};

KerrFluidPoint evaluate_kerr_circular_fluid_point(
    const Metric& metric,
    const Contravariant4& x,
    double mass_M,
    double spin_chi,
    OrbitSense sense);

} // namespace solar::relativity::detail
