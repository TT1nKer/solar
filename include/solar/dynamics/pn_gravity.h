#pragma once

#include "solar/force_model.h"

namespace solar {
namespace dynamics {

// First post-Newtonian (1PN) test-particle acceleration, summed over
// pairs: each body moves in the 1PN field of the others, which is the
// O(1/c^2) expansion of the Schwarzschild geodesic equation:
//
//   a = -(G m / r^2) n [1 - 4 G m / (r c^2) - v^2 / c^2]
//       + 4 (G m / (r^2 c^2)) (v . n) v
//
// with n pointing from the accelerated body toward the source. This is the
// "gradual Newton -> GR" correction layer of the collapse pipeline: it
// vanishes as c -> infinity, reproduces the classical perihelion
// precession (6 pi G M / (a (1-e^2) c^2) per orbit), and is verified
// against Solar's exact Schwarzschild geodesic integrator with O(eps^2)
// residuals (eps = G M / (r c^2)).
//
// Add it alongside a Newtonian force model; the per-particle compactness
// policy (blending windows) belongs to the collapse driver, not to this
// force law.
class PostNewtonianGravity : public ForceModel {
public:
    // speed_of_light_km_s defaults to constants::C_LIGHT; tests override it
    // to verify the c -> infinity reduction to Newtonian gravity.
    explicit PostNewtonianGravity(
        double speed_of_light_km_s = 0.0);

    std::string name() const override { return "gr_1pn_test_particle"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

    // 1PN correction to the potential energy is not tracked here; the
    // Newtonian potential remains the energy diagnostic in the driver.

private:
    double speed_of_light_km_s_ = 0.0;
};

} // namespace dynamics
} // namespace solar
