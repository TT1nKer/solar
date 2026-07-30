#pragma once

#include "solar/relativity/kerr_bl_metric.h"

namespace solar::relativity {

struct KerrConstants {
    double E;
    double Lz;
    double Q;
    double mass_sq;
};

KerrConstants evaluate_kerr_constants(
    const KerrBoyerLindquistMetric& metric,
    const PhaseSpaceState& state,
    GeodesicKind kind);

} // namespace solar::relativity
