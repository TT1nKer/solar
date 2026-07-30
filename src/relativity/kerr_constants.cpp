#include "solar/relativity/kerr_constants.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity {

KerrConstants evaluate_kerr_constants(
    const KerrBoyerLindquistMetric& metric,
    const PhaseSpaceState& state,
    GeodesicKind kind) {
    if (!is_valid_geodesic_kind(kind)) {
        throw std::invalid_argument(
            "Kerr constants geodesic kind is not recognized");
    }
    if (!std::isfinite(state.affine) ||
        !state.x.v.all_finite() ||
        !state.p.v.all_finite()) {
        throw std::domain_error(
            "Kerr constants require a finite phase-space state");
    }
    if (!metric.valid_point(state.x)) {
        throw std::domain_error(
            "Kerr constants state is outside the Kerr BL domain");
    }

    const double energy = -state.p.v[0];
    const double axial_angular_momentum = state.p.v[3];
    const double polar_momentum = state.p.v[2];
    const double theta = state.x.v[2];
    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    const double sin_squared = sin_theta * sin_theta;
    const double cos_squared = cos_theta * cos_theta;
    if (!std::isfinite(sin_squared) || sin_squared == 0.0 ||
        !std::isfinite(cos_squared)) {
        throw std::domain_error(
            "Kerr constants are singular on the BL polar axis");
    }

    const double mass_squared =
        kind == GeodesicKind::Null ? 0.0 : 1.0;
    const double spin = metric.spin_length();
    const double carter =
        polar_momentum * polar_momentum +
        cos_squared *
            (spin * spin *
                 (mass_squared - energy * energy) +
             axial_angular_momentum *
                 axial_angular_momentum /
                 sin_squared);
    if (!std::isfinite(energy) ||
        !std::isfinite(axial_angular_momentum) ||
        !std::isfinite(carter)) {
        throw std::domain_error(
            "Kerr constants evaluation is non-finite");
    }

    return KerrConstants{
        energy,
        axial_angular_momentum,
        carter,
        mass_squared,
    };
}

} // namespace solar::relativity
