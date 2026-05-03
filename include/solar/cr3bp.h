#pragma once

#include "solar/body.h"
#include <vector>
#include <string>

namespace solar {
namespace cr3bp {

// CR3BP system definition
struct System {
    std::string name;
    double mu;          // mass parameter = m2/(m1+m2)
    double L_star;      // characteristic length (km)
    double T_star;      // characteristic time (s)
    double V_star;      // characteristic velocity (km/s) = L_star/T_star
    double mu1, mu2;    // GM of primary and secondary (km^3/s^2)
};

// State in normalized rotating frame
struct CR3BPState {
    double x, y, z;
    double xdot, ydot, zdot;
};

// NOTE: State Transition Matrix (CR3BPStateSTM) was removed in this revision.
// The implementation in src/cr3bp.cpp had a numerical bug — STM values disagreed
// with finite-difference Jacobians by ~100x. The Halo orbit corrector uses
// finite differences instead. If you need STM propagation, please:
//   1. Re-implement and validate against finite differences before use
//   2. See docs/validation/06_halo_orbit_jacobi.md for the bug discussion

struct LagrangePoints {
    CR3BPState L1, L2, L3, L4, L5;
};

struct HaloOrbitParams {
    int lpoint;         // 1 or 2
    double Az;          // z-amplitude (normalized)
    int family;         // 1 = Northern, -1 = Southern
};

struct HaloOrbit {
    std::vector<CR3BPState> trajectory;
    double period;
    double jacobi_constant;
    double stability_index;
    CR3BPState initial_state;
    bool converged;
    int iterations;
};

// Create a system from body names (e.g., "Sun", "Earth")
System make_system(const std::string& primary, const std::string& secondary);

// CR3BP pseudo-potential
double pseudo_potential(double x, double y, double z, double mu);

// CR3BP equations of motion (returns derivatives)
CR3BPState cr3bp_derivatives(const CR3BPState& s, double mu);

// Jacobi constant: C = 2*Omega - v^2
double jacobi_constant(const CR3BPState& s, double mu);

// Lagrange points
LagrangePoints compute_lagrange_points(double mu);

// Propagate CR3BP state
std::vector<CR3BPState> propagate_cr3bp(
    const CR3BPState& ic, double mu, double t_span,
    double dt_output = 0.01, double atol = 1e-12, double rtol = 1e-12);

// Propagate until y=0 crossing (half-period), returns time
double propagate_to_y_crossing(
    CR3BPState& state, double mu, double t_max,
    double atol = 1e-12, double rtol = 1e-12);

// Richardson 3rd-order halo orbit initial conditions
CR3BPState richardson_halo_ic(double mu, int lpoint, double Az, int family = 1);

// Compute halo orbit with differential correction
HaloOrbit compute_halo_orbit(const System& sys, const HaloOrbitParams& params,
                              double tol = 1e-10, int max_iter = 50);

// Unit conversions
State to_physical(const CR3BPState& s, const System& sys);
CR3BPState to_normalized(const State& s, const System& sys);

} // namespace cr3bp
} // namespace solar
