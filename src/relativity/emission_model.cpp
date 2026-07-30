#include "solar/relativity/emission_model.h"

#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

void require_nonnegative_finite(
    double value,
    const char* message) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(message);
    }
}

} // namespace

TransferCoefficients VacuumEmission::coefficients(
    const FluidSample&,
    const Covariant4&,
    double) const {
    return {};
}

GreyEmission::GreyEmission(
    double emissivity_per_density,
    double absorption_per_density)
    : emissivity_per_density_(emissivity_per_density),
      absorption_per_density_(absorption_per_density) {
    require_nonnegative_finite(
        emissivity_per_density_,
        "grey emissivity per density must be finite and non-negative");
    require_nonnegative_finite(
        absorption_per_density_,
        "grey absorption per density must be finite and non-negative");
}

TransferCoefficients GreyEmission::coefficients(
    const FluidSample& fluid,
    const Covariant4& photon_p,
    double observer_frequency) const {
    if (!fluid.valid) {
        return {};
    }
    const double normalized_emitter_frequency =
        -covector_vector_pairing(
            photon_p, fluid.four_velocity);
    const double emitter_frequency =
        observer_frequency * normalized_emitter_frequency;
    const double comoving_emissivity =
        emissivity_per_density_ * fluid.density;
    const double comoving_absorption =
        absorption_per_density_ * fluid.density;
    return TransferCoefficients{
        comoving_emissivity /
            (emitter_frequency * emitter_frequency),
        emitter_frequency * comoving_absorption,
    };
}

DebugPaintEmission::DebugPaintEmission(
    double invariant_emissivity,
    double invariant_absorption)
    : coefficients_{
          invariant_emissivity,
          invariant_absorption} {
    require_nonnegative_finite(
        coefficients_.invariant_emissivity,
        "debug invariant emissivity must be finite and non-negative");
    require_nonnegative_finite(
        coefficients_.invariant_absorption,
        "debug invariant absorption must be finite and non-negative");
}

TransferCoefficients DebugPaintEmission::coefficients(
    const FluidSample& fluid,
    const Covariant4&,
    double) const {
    return fluid.valid
               ? coefficients_
               : TransferCoefficients{};
}

} // namespace solar::relativity
