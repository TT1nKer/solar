#pragma once

#include "solar/relativity/metric.h"

namespace solar::relativity {

struct FluidSample {
    bool valid = false;
    double density = 0.0;
    double temperature = 0.0;
    Contravariant4 four_velocity;
};

class FluidModel {
public:
    virtual ~FluidModel() = default;

    virtual FluidSample sample(
        const Metric& metric,
        const Contravariant4& x) const = 0;
};

class VacuumFluid final : public FluidModel {
public:
    FluidSample sample(
        const Metric& metric,
        const Contravariant4& x) const override;
};

} // namespace solar::relativity
