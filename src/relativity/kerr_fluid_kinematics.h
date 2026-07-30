#pragma once

#include "solar/relativity/fluid_model.h"

namespace solar::relativity::detail {

struct KerrFluidLocation {
    Contravariant4 boyer_lindquist_position;
    double radius;
    double theta;
    double equatorial_height;
};

KerrFluidLocation locate_kerr_fluid_point(
    const Metric& metric,
    const Contravariant4& x,
    double mass_M,
    double spin_chi);

Contravariant4 evaluate_kerr_circular_four_velocity(
    const Metric& metric,
    const Contravariant4& x,
    const KerrFluidLocation& location,
    double mass_M,
    double spin_chi,
    OrbitSense sense);

} // namespace solar::relativity::detail
