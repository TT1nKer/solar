#pragma once

#include "solar/relativity/kerr_orbits.h"
#include "solar/relativity/metric.h"

#include <optional>

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

struct AnalyticCircularDiskConfig {
    double mass_M;
    double spin_chi;
    OrbitSense sense = OrbitSense::Prograde;
    std::optional<double> inner_radius_M;
    double outer_radius_M;
    double density_scale;
    double temperature_scale;
    double density_power = 0.0;
    double surface_height_tolerance = 1.0e-8;
};

/**
 * Equatorial circular disk with a controlled zero-torque temperature shape.
 *
 * Density and temperature use caller-selected model units. This is not a
 * GRMHD or full Page-Thorne disk.
 */
class AnalyticCircularDiskFluid final : public FluidModel {
public:
    explicit AnalyticCircularDiskFluid(
        AnalyticCircularDiskConfig config);

    FluidSample sample(
        const Metric& metric,
        const Contravariant4& x) const override;

    double inner_radius() const noexcept {
        return inner_radius_M_;
    }
    double outer_radius() const noexcept {
        return config_.outer_radius_M;
    }

private:
    AnalyticCircularDiskConfig config_;
    double inner_radius_M_;
};

} // namespace solar::relativity
