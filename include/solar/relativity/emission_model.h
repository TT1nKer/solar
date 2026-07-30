#pragma once

#include "solar/relativity/fluid_model.h"
#include "solar/relativity/radiative_transfer.h"

namespace solar::relativity {

class EmissionModel {
public:
    virtual ~EmissionModel() = default;

    virtual TransferCoefficients coefficients(
        const FluidSample& fluid,
        const Covariant4& photon_p,
        double observer_frequency) const = 0;
};

class VacuumEmission final : public EmissionModel {
public:
    TransferCoefficients coefficients(
        const FluidSample& fluid,
        const Covariant4& photon_p,
        double observer_frequency) const override;
};

class GreyEmission final : public EmissionModel {
public:
    GreyEmission(
        double emissivity_per_density,
        double absorption_per_density);

    TransferCoefficients coefficients(
        const FluidSample& fluid,
        const Covariant4& photon_p,
        double observer_frequency) const override;

private:
    double emissivity_per_density_;
    double absorption_per_density_;
};

/**
 * Diagnostic invariant coefficients for sampling and composition tests.
 *
 * This model is deliberately nonphysical and does not define display color.
 */
class DebugPaintEmission final : public EmissionModel {
public:
    DebugPaintEmission(
        double invariant_emissivity,
        double invariant_absorption);

    TransferCoefficients coefficients(
        const FluidSample& fluid,
        const Covariant4& photon_p,
        double observer_frequency) const override;

private:
    TransferCoefficients coefficients_;
};

} // namespace solar::relativity
