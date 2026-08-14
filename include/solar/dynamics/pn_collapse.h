#pragma once

#include "solar/force_model.h"

#include <vector>

namespace solar {
namespace dynamics {

// Per-particle blending driver for the gradual Newton -> 1PN -> GR
// collapse path (milestone 01, Regime B). One force model supplying the
// full acceleration:
//
//   a_i = g_i + w(eps_i) * fade(eps_i) * [ a_1PN(g_i, Phi_i, v_i) - g_i ]
//
// with the field-based 1PN acceleration
//
//   a_1PN = g [1 + 4 Phi / c^2 - v^2 / c^2] - 4 (g . v) v / c^2
//
// and per-particle compactness eps_i = -Phi_i / c^2 = G M_enc / (r_i c^2)
// from the same Barnes-Hut monopole walk. The blending weight ramps the
// correction on smoothly across the window (smoothstep):
//
//   eps <= eps_lo          -> pure Newtonian (w = 0)
//   eps_lo < eps < eps_hi  -> Newton + 1PN, smoothstep-weighted
//   eps >= eps_hi          -> full 1PN (w = 1), subject to fade()
//
// fade() (same policy as PostNewtonianGravity) suppresses the correction
// beyond the 1PN validity boundary fade_hi so dense phases never feed
// unphysical velocity-coupled terms to the integrator. The eps >= eps_hi
// hand-off to the LTB compact stage belongs to the next milestone and is
// exposed here through the handoff_candidates() diagnostic.
//
// Degenerate windows: eps_lo = 0, eps_hi = 0 selects the full 1PN force
// everywhere (eps > 0 -> w = 1), matching an un-blended
// PostNewtonianGravity with the same config; a very large eps_lo selects
// pure Newtonian tree gravity. Both reductions are verified in
// test_collapse_blend.
class PnCollapseForce : public ForceModel {
public:
    struct Config {
        double opening_angle = 0.7;   // theta for the monopole walk
        double softening_km = 0.0;    // Plummer softening length (km)
        double eps_lo = 1.0e-4;       // below: pure Newtonian
        double eps_hi = 0.05;         // above: full 1PN (before fade)
        double fade_hi = 0.1;         // validity fade threshold; 0 = off
    };

    PnCollapseForce() = default;
    explicit PnCollapseForce(Config config);
    PnCollapseForce(Config config, double speed_of_light_km_s);

    std::string name() const override { return "pn_collapse_blend"; }

    void compute(
        const std::vector<Body>& bodies,
        double time,
        std::vector<Vec3>& acc) const override;

    // Exact softened Newtonian potential (same direct-sum form as
    // BarnesHutGravity::potential_energy) so the driver's energy
    // diagnostics stay unpolluted by the tree or the PN layer.
    double potential_energy(
        const std::vector<Body>& bodies,
        double time) const override;

private:
    Config config_;
    double speed_of_light_km_s_ = 0.0;
};

// Diagnostic: exact per-particle compactness eps_i = -Phi_i / c^2 via
// direct O(N^2) softened summation (no tree truncation).
std::vector<double> per_particle_compactness(
    const std::vector<Body>& bodies,
    double softening_km,
    double speed_of_light_km_s);

// Diagnostic: indices of particles with eps_i >= eps_hi — the terminal
// core to hand off to the LTB compact stage.
std::vector<int> handoff_candidates(
    const std::vector<Body>& bodies,
    double softening_km,
    double speed_of_light_km_s,
    double eps_hi);

} // namespace dynamics
} // namespace solar
