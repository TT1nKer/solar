#pragma once

#include "solar/relativity/kerr_schild_metric.h"

namespace solar::relativity {

/**
 * Complete canonical transform between safe exterior Boyer-Lindquist and
 * ingoing Cartesian Kerr-Schild coordinates.
 *
 * The transform is deliberately restricted to a finite exterior overlap.
 * Horizon crossing uses states already expressed in Kerr-Schild coordinates.
 * Covariant canonical momenta are transformed with the inverse-transpose
 * Jacobian; they are never treated as spatial vectors.
 */
class KerrChartTransform {
public:
    KerrChartTransform(
        double mass_M,
        double spin_chi,
        double overlap_margin_fraction = 1.0e-4);

    Contravariant4 position_to_kerr_schild(
        const Contravariant4& boyer_lindquist) const;
    Contravariant4 position_to_boyer_lindquist(
        const Contravariant4& kerr_schild) const;

    Mat4 boyer_lindquist_to_kerr_schild_jacobian(
        const Contravariant4& boyer_lindquist) const;
    Mat4 kerr_schild_to_boyer_lindquist_jacobian(
        const Contravariant4& kerr_schild) const;

    PhaseSpaceState state_to_kerr_schild(
        const PhaseSpaceState& boyer_lindquist) const;
    PhaseSpaceState state_to_boyer_lindquist(
        const PhaseSpaceState& kerr_schild) const;

private:
    KerrSchildCartesianMetric kerr_schild_metric_;
    double overlap_margin_M_;
    double minimum_overlap_radius_M_;
};

} // namespace solar::relativity
