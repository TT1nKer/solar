#pragma once

#include "solar/relativity/metric.h"

namespace solar::relativity {

/**
 * Kerr spacetime in ingoing Cartesian Kerr-Schild coordinates (t,x,y,z).
 *
 * The chart is regular on the outer horizon and supports positive implicit
 * radius down to the ring/zero-radius boundary. Parameters use geometrized
 * units with mass_M > 0 and dimensionless spin abs(spin_chi) < 1.
 *
 * Validation:
 *   tests/relativity/test_kerr_schild.cpp
 *   tests/relativity/test_kerr_schild_derivatives.cpp
 */
class KerrSchildCartesianMetric final : public Metric {
public:
    KerrSchildCartesianMetric(
        double mass_M,
        double spin_chi);

    Chart chart() const noexcept override {
        return Chart::KerrSchildCartesian;
    }
    std::string name() const override {
        return "kerr-schild-cartesian";
    }

    Mat4 covariant(const Contravariant4& x) const override;
    Mat4 contravariant(const Contravariant4& x) const override;
    std::array<Mat4, 4>
    contravariant_derivatives(
        const Contravariant4& x) const override;
    bool valid_point(
        const Contravariant4& x) const noexcept override;

    double radial_coordinate(
        const Contravariant4& x) const;
    Vec3 radial_coordinate_gradient(
        const Contravariant4& x) const;

    double mass() const noexcept { return mass_M_; }
    double spin_chi() const noexcept { return spin_chi_; }
    double spin_length() const noexcept { return spin_a_M_; }
    double outer_horizon_radius() const noexcept {
        return outer_horizon_M_;
    }
    double inner_horizon_radius() const noexcept {
        return inner_horizon_M_;
    }

private:
    double mass_M_;
    double spin_chi_;
    double spin_a_M_;
    double outer_horizon_M_;
    double inner_horizon_M_;
    double minimum_radius_M_;
};

double kerr_schild_stationary_energy(
    const PhaseSpaceState& state);

double kerr_schild_axial_angular_momentum(
    const PhaseSpaceState& state);

} // namespace solar::relativity
