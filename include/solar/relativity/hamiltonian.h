#pragma once

#include "solar/numerics/dopri5.h"
#include "solar/relativity/metric.h"

namespace solar::relativity {

struct PhaseSpaceDerivative {
    Contravariant4 dx;
    Covariant4 dp;
};

double hamiltonian(
    const Metric& metric,
    const PhaseSpaceState& state);

double hamiltonian_constraint_error(
    const Metric& metric,
    const PhaseSpaceState& state,
    GeodesicKind kind);

class HamiltonGeodesicRhs {
public:
    explicit HamiltonGeodesicRhs(const Metric& metric) noexcept;

    PhaseSpaceDerivative operator()(
        const PhaseSpaceState& state) const;

private:
    const Metric* metric_;
};

numerics::StateN<8> pack_phase_space(
    const PhaseSpaceState& state);

PhaseSpaceState unpack_phase_space(
    double affine,
    const numerics::StateN<8>& packed);

} // namespace solar::relativity
