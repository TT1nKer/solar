#pragma once

#include "solar/force_model.h"

namespace solar {
namespace dynamics {

// First post-Newtonian (1PN) correction in the standard PN gauge,
// expressed through the Newtonian field g and potential Phi of the
// other bodies:
//
//   a = g [1 + 4 Phi / c^2 - v^2 / c^2] - 4 (g . v) v / c^2
//
// For a single point source this is exactly the classical test-particle
// acceleration
//
//   a = -(G m / r^2) n [1 - 4 G m / (r c^2) - v^2 / c^2]
//       + 4 (G m / (r^2 c^2)) (v . n) v,
//
// so the analytic anchors carry over unchanged: the force reduces to
// Newtonian gravity as c -> infinity, matches the radial spot check, and
// reproduces the perihelion precession 6 pi G M / (a (1-e^2) c^2) per
// orbit. For an extended cloud the field form additionally reproduces
// the mean-field 1PN dynamics (the surface shell of a spherical cloud
// follows r'' = -(G M / r^2) [1 - 4 G M / (r c^2) - 5 v^2 / c^2], with
// the enclosed mass M in both corrections), which the pairwise
// test-particle form cannot: its static term scales per-pair and misses
// the enclosed-mass O(G M / (r c^2)) term. The spherical-collapse
// regression anchors on this mean-field limit.
//
// Field and potential are evaluated with the same Barnes-Hut octree as
// BarnesHutGravity. Config::opening_angle == 0 disables the multipole
// approximation (exact O(N^2) summation, matching the analytic anchors
// bit for bit); pass the same theta/softening as the companion Newtonian
// model for large-N collapse runs. The model contributes the full
// acceleration (base plus correction), so it replaces rather than
// supplements a Newtonian force model.
//
// Per-particle fade: the 1PN expansion is only valid for small local
// compactness eps = -Phi / c^2 = G M_enc / (r c^2). The contribution is
// scaled by 1 - smoothstep(eps / fade_hi) for eps < fade_hi and zero
// beyond, so dense phases never feed unphysical velocity-coupled terms
// into the integrator. The default (0) keeps the pure standard-gauge
// force so the analytic anchors hold exactly; the collapse blending
// driver opts into a non-zero threshold where it needs it.
class PostNewtonianGravity : public ForceModel {
public:
    struct Config {
        double opening_angle = 0.0;   // theta; 0 = exact direct sum
        double softening_km = 0.0;    // Plummer softening length (km)
        double fade_hi = 0.0;         // compactness fade threshold; 0 = off (default)
    };

    // speed_of_light_km_s defaults to constants::C_LIGHT; tests override
    // it to verify the c -> infinity reduction to Newtonian gravity.
    explicit PostNewtonianGravity(
        double speed_of_light_km_s = 0.0);
    PostNewtonianGravity(Config config, double speed_of_light_km_s = 0.0);

    std::string name() const override { return "gr_1pn_field"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

    // The 1PN correction to the potential energy is not tracked here;
    // the Newtonian potential remains the energy diagnostic in the
    // driver.

private:
    Config config_;
    double speed_of_light_km_s_ = 0.0;
};

} // namespace dynamics
} // namespace solar