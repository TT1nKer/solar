#include "solar/integrator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace solar {

static bool is_finite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

static void require_finite_accelerations(
    const std::vector<Vec3>& accelerations, size_t expected_size, const char* integrator_name) {
    if (accelerations.size() != expected_size)
        throw std::runtime_error(std::string(integrator_name) + ": acceleration size mismatch");
    for (const auto& acceleration : accelerations) {
        if (!is_finite(acceleration))
            throw std::runtime_error(std::string(integrator_name) + ": non-finite acceleration");
    }
}

std::vector<State> rk4_step(
    const std::vector<State>& states,
    const std::vector<double>& masses,
    const AccelFunc& accel,
    double time,
    double dt)
{
    size_t n = states.size();
    if (masses.size() != n) throw std::invalid_argument("rk4_step: state/mass size mismatch");
    if (!std::isfinite(time) || !std::isfinite(dt) || dt == 0.0)
        throw std::invalid_argument("rk4_step: time and non-zero dt must be finite");

    // k1
    auto a1 = accel(time, states, masses);
    require_finite_accelerations(a1, n, "rk4_step");
    std::vector<State> s1(n);
    for (size_t i = 0; i < n; ++i) {
        s1[i].pos = states[i].pos + states[i].vel * (dt * 0.5);
        s1[i].vel = states[i].vel + a1[i] * (dt * 0.5);
    }

    // k2
    auto a2 = accel(time + 0.5 * dt, s1, masses);
    require_finite_accelerations(a2, n, "rk4_step");
    std::vector<State> s2(n);
    for (size_t i = 0; i < n; ++i) {
        s2[i].pos = states[i].pos + s1[i].vel * (dt * 0.5);
        s2[i].vel = states[i].vel + a2[i] * (dt * 0.5);
    }

    // k3
    auto a3 = accel(time + 0.5 * dt, s2, masses);
    require_finite_accelerations(a3, n, "rk4_step");
    std::vector<State> s3(n);
    for (size_t i = 0; i < n; ++i) {
        s3[i].pos = states[i].pos + s2[i].vel * dt;
        s3[i].vel = states[i].vel + a3[i] * dt;
    }

    // k4
    auto a4 = accel(time + dt, s3, masses);
    require_finite_accelerations(a4, n, "rk4_step");

    // Combine: state += (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    std::vector<State> result(n);
    for (size_t i = 0; i < n; ++i) {
        result[i].pos = states[i].pos +
            (states[i].vel + 2.0 * s1[i].vel + 2.0 * s2[i].vel + s3[i].vel) * (dt / 6.0);
        result[i].vel = states[i].vel +
            (a1[i] + 2.0 * a2[i] + 2.0 * a3[i] + a4[i]) * (dt / 6.0);
    }

    return result;
}

void verlet_step(
    std::vector<State>& states,
    const std::vector<double>& masses,
    const AccelFunc& accel,
    std::vector<Vec3>& prev_accels,
    double time,
    double dt)
{
    size_t n = states.size();
    if (masses.size() != n || prev_accels.size() != n)
        throw std::invalid_argument("verlet_step: state, mass, and acceleration size mismatch");
    if (!std::isfinite(time) || !std::isfinite(dt) || dt == 0.0)
        throw std::invalid_argument("verlet_step: time and non-zero dt must be finite");

    // Update positions: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    for (size_t i = 0; i < n; ++i) {
        states[i].pos += states[i].vel * dt + prev_accels[i] * (0.5 * dt * dt);
    }

    // Compute new accelerations at updated positions
    auto new_accels = accel(time + dt, states, masses);
    require_finite_accelerations(new_accels, n, "verlet_step");

    // Update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    for (size_t i = 0; i < n; ++i) {
        states[i].vel += (prev_accels[i] + new_accels[i]) * (0.5 * dt);
    }

    prev_accels = new_accels;
}

// ============================================================
// Dormand-Prince RK5(4) — DOPRI5
// Reference: Dormand & Prince (1980), J. Comp. Appl. Math. 6
// Hairer, Norsett, Wanner "Solving ODEs I", Table 5.2
// ============================================================

// Butcher tableau coefficients
static constexpr double dp_a21 = 1.0/5.0;
static constexpr double dp_a31 = 3.0/40.0,       dp_a32 = 9.0/40.0;
static constexpr double dp_a41 = 44.0/45.0,      dp_a42 = -56.0/15.0,     dp_a43 = 32.0/9.0;
static constexpr double dp_a51 = 19372.0/6561.0,  dp_a52 = -25360.0/2187.0, dp_a53 = 64448.0/6561.0, dp_a54 = -212.0/729.0;
static constexpr double dp_a61 = 9017.0/3168.0,   dp_a62 = -355.0/33.0,    dp_a63 = 46732.0/5247.0,  dp_a64 = 49.0/176.0,    dp_a65 = -5103.0/18656.0;

// 5th-order weights (b)
static constexpr double dp_b1 = 35.0/384.0,  dp_b3 = 500.0/1113.0, dp_b4 = 125.0/192.0;
static constexpr double dp_b5 = -2187.0/6784.0, dp_b6 = 11.0/84.0;
// b7 = 0 (FSAL: k7 = f(t+h, y_new))

// 4th-order weights (b*) for error estimation
static constexpr double dp_e1 = 71.0/57600.0,    dp_e3 = -71.0/16695.0, dp_e4 = 71.0/1920.0;
static constexpr double dp_e5 = -17253.0/339200.0, dp_e6 = 22.0/525.0,   dp_e7 = -1.0/40.0;
// e_i = b_i - b*_i (error coefficients)

// c nodes
static constexpr double dp_c2 = 1.0/5.0, dp_c3 = 3.0/10.0, dp_c4 = 4.0/5.0;
static constexpr double dp_c5 = 8.0/9.0, dp_c6 = 1.0;

AdaptiveResult dopri5_step(
    const std::vector<State>& states,
    const std::vector<double>& masses,
    const AccelFunc& accel,
    double time,
    double dt,
    double atol,
    double rtol)
{
    size_t n = states.size();
    if (masses.size() != n) throw std::invalid_argument("dopri5_step: state/mass size mismatch");
    if (!std::isfinite(time) || !std::isfinite(dt) || dt == 0.0 ||
        !std::isfinite(atol) || !std::isfinite(rtol) || atol <= 0.0 || rtol < 0.0)
        throw std::invalid_argument("dopri5_step: invalid time, step, or tolerance");

    // Helper: compute state at s0 + dt*(sum of a_ij * k_j) for velocity,
    // and s0.pos + dt*(sum of a_ij * v_j) for position
    // k stages store {velocity, acceleration}
    struct Stage { Vec3 vel; Vec3 acc; };

    // k1: evaluate at t, y
    auto acc1 = accel(time, states, masses);
    require_finite_accelerations(acc1, n, "dopri5_step");
    std::vector<Stage> k1(n);
    for (size_t i = 0; i < n; ++i) {
        k1[i] = {states[i].vel, acc1[i]};
    }

    // k2: evaluate at t + c2*dt
    std::vector<State> s2(n);
    for (size_t i = 0; i < n; ++i) {
        s2[i].pos = states[i].pos + k1[i].vel * (dp_a21 * dt);
        s2[i].vel = states[i].vel + k1[i].acc * (dp_a21 * dt);
    }
    auto acc2 = accel(time + dp_c2 * dt, s2, masses);
    require_finite_accelerations(acc2, n, "dopri5_step");
    std::vector<Stage> k2(n);
    for (size_t i = 0; i < n; ++i) k2[i] = {s2[i].vel, acc2[i]};

    // k3
    std::vector<State> s3(n);
    for (size_t i = 0; i < n; ++i) {
        s3[i].pos = states[i].pos + (k1[i].vel * dp_a31 + k2[i].vel * dp_a32) * dt;
        s3[i].vel = states[i].vel + (k1[i].acc * dp_a31 + k2[i].acc * dp_a32) * dt;
    }
    auto acc3 = accel(time + dp_c3 * dt, s3, masses);
    require_finite_accelerations(acc3, n, "dopri5_step");
    std::vector<Stage> k3(n);
    for (size_t i = 0; i < n; ++i) k3[i] = {s3[i].vel, acc3[i]};

    // k4
    std::vector<State> s4(n);
    for (size_t i = 0; i < n; ++i) {
        s4[i].pos = states[i].pos + (k1[i].vel * dp_a41 + k2[i].vel * dp_a42 + k3[i].vel * dp_a43) * dt;
        s4[i].vel = states[i].vel + (k1[i].acc * dp_a41 + k2[i].acc * dp_a42 + k3[i].acc * dp_a43) * dt;
    }
    auto acc4 = accel(time + dp_c4 * dt, s4, masses);
    require_finite_accelerations(acc4, n, "dopri5_step");
    std::vector<Stage> k4(n);
    for (size_t i = 0; i < n; ++i) k4[i] = {s4[i].vel, acc4[i]};

    // k5
    std::vector<State> s5(n);
    for (size_t i = 0; i < n; ++i) {
        s5[i].pos = states[i].pos + (k1[i].vel * dp_a51 + k2[i].vel * dp_a52
                    + k3[i].vel * dp_a53 + k4[i].vel * dp_a54) * dt;
        s5[i].vel = states[i].vel + (k1[i].acc * dp_a51 + k2[i].acc * dp_a52
                    + k3[i].acc * dp_a53 + k4[i].acc * dp_a54) * dt;
    }
    auto acc5 = accel(time + dp_c5 * dt, s5, masses);
    require_finite_accelerations(acc5, n, "dopri5_step");
    std::vector<Stage> k5(n);
    for (size_t i = 0; i < n; ++i) k5[i] = {s5[i].vel, acc5[i]};

    // k6
    std::vector<State> s6(n);
    for (size_t i = 0; i < n; ++i) {
        s6[i].pos = states[i].pos + (k1[i].vel * dp_a61 + k2[i].vel * dp_a62
                    + k3[i].vel * dp_a63 + k4[i].vel * dp_a64 + k5[i].vel * dp_a65) * dt;
        s6[i].vel = states[i].vel + (k1[i].acc * dp_a61 + k2[i].acc * dp_a62
                    + k3[i].acc * dp_a63 + k4[i].acc * dp_a64 + k5[i].acc * dp_a65) * dt;
    }
    auto acc6 = accel(time + dp_c6 * dt, s6, masses);
    require_finite_accelerations(acc6, n, "dopri5_step");
    std::vector<Stage> k6(n);
    for (size_t i = 0; i < n; ++i) k6[i] = {s6[i].vel, acc6[i]};

    // 5th-order solution (y_{n+1})
    std::vector<State> y_new(n);
    for (size_t i = 0; i < n; ++i) {
        y_new[i].pos = states[i].pos + (k1[i].vel * dp_b1 + k3[i].vel * dp_b3
                       + k4[i].vel * dp_b4 + k5[i].vel * dp_b5 + k6[i].vel * dp_b6) * dt;
        y_new[i].vel = states[i].vel + (k1[i].acc * dp_b1 + k3[i].acc * dp_b3
                       + k4[i].acc * dp_b4 + k5[i].acc * dp_b5 + k6[i].acc * dp_b6) * dt;
    }

    // k7 (FSAL: evaluate at new state for error estimate)
    auto acc7 = accel(time + dt, y_new, masses);
    require_finite_accelerations(acc7, n, "dopri5_step");

    // Error estimate: difference between 5th and 4th order solutions
    // e = dt * sum(e_i * k_i)  where e_i = b_i - b*_i
    double err_max = 0.0;
    for (size_t i = 0; i < n; ++i) {
        Vec3 err_pos = (k1[i].vel * dp_e1 + k3[i].vel * dp_e3 + k4[i].vel * dp_e4
                       + k5[i].vel * dp_e5 + k6[i].vel * dp_e6 + Vec3{y_new[i].vel.x, y_new[i].vel.y, y_new[i].vel.z} * dp_e7) * dt;
        // Use k7 velocity for the error e7 term (FSAL)
        err_pos = (k1[i].vel * dp_e1 + k3[i].vel * dp_e3 + k4[i].vel * dp_e4
                  + k5[i].vel * dp_e5 + k6[i].vel * dp_e6) * dt;
        // k7 contribution to error for position
        Vec3 k7_vel = {y_new[i].vel.x, y_new[i].vel.y, y_new[i].vel.z};
        err_pos += k7_vel * (dp_e7 * dt);

        Vec3 err_vel = (k1[i].acc * dp_e1 + k3[i].acc * dp_e3 + k4[i].acc * dp_e4
                       + k5[i].acc * dp_e5 + k6[i].acc * dp_e6 + acc7[i] * dp_e7) * dt;

        // Scale by tolerance: err_scaled = err / (atol + rtol * max(|y|, |y_new|))
        double sc_pos = atol + rtol * std::max(states[i].pos.norm(), y_new[i].pos.norm());
        double sc_vel = atol + rtol * std::max(states[i].vel.norm(), y_new[i].vel.norm());

        double ep = err_pos.norm() / sc_pos;
        double ev = err_vel.norm() / sc_vel;
        if (!std::isfinite(ep) || !std::isfinite(ev))
            throw std::runtime_error("dopri5_step: non-finite state or error estimate");
        err_max = std::max(err_max, std::max(ep, ev));
    }

    // Step size control
    // err_max should be < 1 for acceptance
    AdaptiveResult result;
    result.states = std::move(y_new);
    result.dt_used = dt;
    result.error = err_max;
    result.accepted = (err_max <= 1.0);

    // Recommended next step: dt_new = dt * safety * (1/err)^(1/5)
    double safety = 0.9;
    double order = 5.0;
    double factor;
    if (err_max > 1e-30) {
        factor = safety * std::pow(1.0 / err_max, 1.0 / order);
    } else {
        factor = 5.0; // error essentially zero, can increase a lot
    }
    factor = std::max(0.2, std::min(5.0, factor)); // clamp growth/shrink
    result.dt_next = dt * factor;

    return result;
}

// ============================================================
// Generic DOPRI5 for arbitrary ODE
// Reuses the same Butcher tableau as above
// ============================================================

GenericAdaptiveResult dopri5_generic_step(
    const std::vector<double>& y, double t, double dt,
    const GenericODEFunc& f,
    double atol, double rtol)
{
    size_t n = y.size();

    if (n == 0) throw std::invalid_argument("dopri5_generic_step: state must not be empty");
    if (!std::isfinite(t) || !std::isfinite(dt) || dt == 0.0 ||
        !std::isfinite(atol) || !std::isfinite(rtol) || atol <= 0.0 || rtol < 0.0)
        throw std::invalid_argument("dopri5_generic_step: invalid time, step, or tolerance");
    for (double value : y) {
        if (!std::isfinite(value))
            throw std::invalid_argument("dopri5_generic_step: state must be finite");
    }

    auto evaluate = [&](double stage_time, const std::vector<double>& stage_state) {
        auto derivative = f(stage_time, stage_state);
        if (derivative.size() != n)
            throw std::runtime_error("dopri5_generic_step: derivative size mismatch");
        for (double value : derivative) {
            if (!std::isfinite(value))
                throw std::runtime_error("dopri5_generic_step: non-finite derivative");
        }
        return derivative;
    };

    auto k1 = evaluate(t, y);

    std::vector<double> ys(n);

    // k2
    for (size_t i = 0; i < n; ++i) ys[i] = y[i] + dp_a21 * dt * k1[i];
    auto k2 = evaluate(t + dp_c2 * dt, ys);

    // k3
    for (size_t i = 0; i < n; ++i) ys[i] = y[i] + dt * (dp_a31*k1[i] + dp_a32*k2[i]);
    auto k3 = evaluate(t + dp_c3 * dt, ys);

    // k4
    for (size_t i = 0; i < n; ++i) ys[i] = y[i] + dt * (dp_a41*k1[i] + dp_a42*k2[i] + dp_a43*k3[i]);
    auto k4 = evaluate(t + dp_c4 * dt, ys);

    // k5
    for (size_t i = 0; i < n; ++i) ys[i] = y[i] + dt * (dp_a51*k1[i] + dp_a52*k2[i] + dp_a53*k3[i] + dp_a54*k4[i]);
    auto k5 = evaluate(t + dp_c5 * dt, ys);

    // k6
    for (size_t i = 0; i < n; ++i) ys[i] = y[i] + dt * (dp_a61*k1[i] + dp_a62*k2[i] + dp_a63*k3[i] + dp_a64*k4[i] + dp_a65*k5[i]);
    auto k6 = evaluate(t + dp_c6 * dt, ys);

    // 5th-order solution
    std::vector<double> y_new(n);
    for (size_t i = 0; i < n; ++i) {
        y_new[i] = y[i] + dt * (dp_b1*k1[i] + dp_b3*k3[i] + dp_b4*k4[i] + dp_b5*k5[i] + dp_b6*k6[i]);
    }

    // k7 (FSAL)
    auto k7 = evaluate(t + dt, y_new);

    // Error estimate
    double err_max = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double err_i = dt * (dp_e1*k1[i] + dp_e3*k3[i] + dp_e4*k4[i]
                           + dp_e5*k5[i] + dp_e6*k6[i] + dp_e7*k7[i]);
        double sc = atol + rtol * std::max(std::fabs(y[i]), std::fabs(y_new[i]));
        double scaled_error = std::fabs(err_i) / sc;
        if (!std::isfinite(scaled_error))
            throw std::runtime_error("dopri5_generic_step: non-finite state or error estimate");
        err_max = std::max(err_max, scaled_error);
    }

    GenericAdaptiveResult result;
    result.y = std::move(y_new);
    result.dt_used = dt;
    result.error = err_max;
    result.accepted = (err_max <= 1.0);

    double factor;
    if (err_max > 1e-30) {
        factor = 0.9 * std::pow(1.0 / err_max, 0.2);
    } else {
        factor = 5.0;
    }
    factor = std::max(0.2, std::min(5.0, factor));
    result.dt_next = dt * factor;

    return result;
}

} // namespace solar
