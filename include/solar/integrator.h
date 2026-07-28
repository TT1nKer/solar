#pragma once

#include "solar/body.h"
#include <vector>
#include <functional>

namespace solar {

// Acceleration function evaluated at one integrator stage.
// The complete state and stage time are required for velocity- and time-dependent forces.
using AccelFunc = std::function<std::vector<Vec3>(
    double time,
    const std::vector<State>& states,
    const std::vector<double>& masses)>;

// 4th-order Runge-Kutta step for N bodies
// Returns updated states after one timestep
std::vector<State> rk4_step(
    const std::vector<State>& states,
    const std::vector<double>& masses,
    const AccelFunc& accel,
    double time,
    double dt);

// Velocity Verlet (symplectic) step — modifies states in place
// prev_accels: accelerations from previous step (updated on return)
void verlet_step(
    std::vector<State>& states,
    const std::vector<double>& masses,
    const AccelFunc& accel,
    std::vector<Vec3>& prev_accels,
    double time,
    double dt);

// Adaptive Dormand-Prince step result
struct AdaptiveResult {
    std::vector<State> states;  // updated states (higher-order solution)
    double dt_used;             // timestep actually taken
    double dt_next;             // recommended next timestep
    double error;               // estimated relative error norm
    bool accepted;              // true if error < tolerance
};

// Dormand-Prince RK5(4) adaptive step (DOPRI5)
// 7-stage embedded pair: 5th-order solution advances state,
// 4th-order provides error estimate.
// atol: absolute tolerance (km for position, km/s for velocity)
// rtol: relative tolerance
AdaptiveResult dopri5_step(
    const std::vector<State>& states,
    const std::vector<double>& masses,
    const AccelFunc& accel,
    double time,
    double dt,
    double atol = 1e-10,
    double rtol = 1e-10);

// Generic DOPRI5 for arbitrary ODE: dy/dt = f(t, y)
// Works on a flat vector of doubles — not limited to N-body
struct GenericAdaptiveResult {
    std::vector<double> y;
    double dt_used;
    double dt_next;
    double error;
    bool accepted;
};

using GenericODEFunc = std::function<std::vector<double>(double t, const std::vector<double>& y)>;

GenericAdaptiveResult dopri5_generic_step(
    const std::vector<double>& y, double t, double dt,
    const GenericODEFunc& f,
    double atol = 1e-12, double rtol = 1e-12);

} // namespace solar
