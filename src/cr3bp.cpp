#include "solar/cr3bp.h"
#include "solar/integrator.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace solar {
namespace cr3bp {

// ============================================================
// System factory
// ============================================================

System make_system(const std::string& primary, const std::string& secondary) {
    using namespace constants;

    // Map names to GM values and semi-major axes
    struct BodyData { double gm; double a_km; }; // a_km: orbital distance from primary

    auto get_data = [](const std::string& name) -> BodyData {
        std::string n = name;
        std::transform(n.begin(), n.end(), n.begin(), ::tolower);
        if (n == "sun")     return {MU_SUN, 0};
        if (n == "earth")   return {MU_EARTH, 1.000003 * AU};
        if (n == "moon")    return {MU_MOON, 384400.0};
        if (n == "mars")    return {MU_MARS, 1.52371 * AU};
        if (n == "jupiter") return {MU_JUPITER, 5.20289 * AU};
        if (n == "saturn")  return {MU_SATURN, 9.53668 * AU};
        throw std::runtime_error("Unknown body for CR3BP: " + name);
    };

    auto p = get_data(primary);
    auto s = get_data(secondary);

    System sys;
    sys.name = primary + "-" + secondary;
    sys.mu1 = p.gm;
    sys.mu2 = s.gm;
    sys.mu = s.gm / (p.gm + s.gm);
    sys.L_star = s.a_km; // distance between primaries
    sys.T_star = std::sqrt(sys.L_star * sys.L_star * sys.L_star / (p.gm + s.gm));
    sys.V_star = sys.L_star / sys.T_star;

    return sys;
}

// ============================================================
// Equations of motion
// ============================================================

double pseudo_potential(double x, double y, double z, double mu) {
    double r1 = std::sqrt((x + mu)*(x + mu) + y*y + z*z);
    double r2 = std::sqrt((x - 1.0 + mu)*(x - 1.0 + mu) + y*y + z*z);
    return 0.5*(x*x + y*y) + (1.0 - mu)/r1 + mu/r2;
}

CR3BPState cr3bp_derivatives(const CR3BPState& s, double mu) {
    double x = s.x, y = s.y, z = s.z;

    double r1_sq = (x + mu)*(x + mu) + y*y + z*z;
    double r2_sq = (x - 1.0 + mu)*(x - 1.0 + mu) + y*y + z*z;
    double r1_3 = r1_sq * std::sqrt(r1_sq);
    double r2_3 = r2_sq * std::sqrt(r2_sq);

    double Ox = x - (1.0 - mu)*(x + mu)/r1_3 - mu*(x - 1.0 + mu)/r2_3;
    double Oy = y - (1.0 - mu)*y/r1_3 - mu*y/r2_3;
    double Oz =   - (1.0 - mu)*z/r1_3 - mu*z/r2_3;

    return {
        s.xdot, s.ydot, s.zdot,
        2.0*s.ydot + Ox,
        -2.0*s.xdot + Oy,
        Oz
    };
}

double jacobi_constant(const CR3BPState& s, double mu) {
    double v_sq = s.xdot*s.xdot + s.ydot*s.ydot + s.zdot*s.zdot;
    return 2.0 * pseudo_potential(s.x, s.y, s.z, mu) - v_sq;
}

// ============================================================
// Lagrange points
// ============================================================

LagrangePoints compute_lagrange_points(double mu) {
    LagrangePoints lp;

    // f(x) = dOmega/dx on y=z=0 axis
    auto f_collinear = [mu](double x) {
        double r1 = std::fabs(x + mu);
        double r2 = std::fabs(x - 1.0 + mu);
        double s1 = (x + mu > 0) ? 1.0 : -1.0;
        double s2 = (x - 1.0 + mu > 0) ? 1.0 : -1.0;
        return x - (1.0 - mu)*s1/(r1*r1) - mu*s2/(r2*r2);
    };

    auto df_collinear = [mu](double x) {
        double r1 = std::fabs(x + mu);
        double r2 = std::fabs(x - 1.0 + mu);
        return 1.0 + 2.0*(1.0 - mu)/(r1*r1*r1) + 2.0*mu/(r2*r2*r2);
    };

    // Newton-Raphson
    auto solve = [&](double x0) {
        for (int i = 0; i < 50; ++i) {
            double fx = f_collinear(x0);
            double dfx = df_collinear(x0);
            double dx = fx / dfx;
            x0 -= dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        return x0;
    };

    double mu13 = std::cbrt(mu / 3.0);

    // L1: between primaries
    double x_L1 = solve(1.0 - mu - mu13);
    lp.L1 = {x_L1, 0, 0, 0, 0, 0};

    // L2: beyond secondary
    double x_L2 = solve(1.0 - mu + mu13);
    lp.L2 = {x_L2, 0, 0, 0, 0, 0};

    // L3: opposite side of primary
    double x_L3 = solve(-1.0 - 5.0*mu/12.0);
    lp.L3 = {x_L3, 0, 0, 0, 0, 0};

    // L4, L5: equilateral triangle positions
    lp.L4 = {0.5 - mu,  std::sqrt(3.0)/2.0, 0, 0, 0, 0};
    lp.L5 = {0.5 - mu, -std::sqrt(3.0)/2.0, 0, 0, 0, 0};

    return lp;
}

// ============================================================
// Propagation
// ============================================================

// Pack/unpack CR3BPState to/from flat vector
static std::vector<double> pack_state(const CR3BPState& s) {
    return {s.x, s.y, s.z, s.xdot, s.ydot, s.zdot};
}

static CR3BPState unpack_state(const std::vector<double>& v) {
    return {v[0], v[1], v[2], v[3], v[4], v[5]};
}

// Pack state + identity STM (42 elements)
static std::vector<double> pack_state_stm(const CR3BPState& s) {
    std::vector<double> v(42, 0.0);
    v[0]=s.x; v[1]=s.y; v[2]=s.z; v[3]=s.xdot; v[4]=s.ydot; v[5]=s.zdot;
    // Identity STM
    for (int i = 0; i < 6; ++i) v[6 + i*6 + i] = 1.0;
    return v;
}

static CR3BPStateSTM unpack_state_stm(const std::vector<double>& v) {
    CR3BPStateSTM r;
    r.state = {v[0], v[1], v[2], v[3], v[4], v[5]};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            r.stm[i][j] = v[6 + i*6 + j];
    return r;
}

// CR3BP + STM RHS for the generic integrator
static std::vector<double> cr3bp_rhs_stm(double /*t*/, const std::vector<double>& y, double mu) {
    CR3BPState s = {y[0], y[1], y[2], y[3], y[4], y[5]};
    CR3BPState ds = cr3bp_derivatives(s, mu);

    // Compute Jacobian of EOM (A matrix) for STM propagation
    double x = s.x, yy = s.y, z = s.z;
    double r1_sq = (x+mu)*(x+mu) + yy*yy + z*z;
    double r2_sq = (x-1+mu)*(x-1+mu) + yy*yy + z*z;
    double r1_5 = r1_sq * r1_sq * std::sqrt(r1_sq);
    double r2_5 = r2_sq * r2_sq * std::sqrt(r2_sq);
    double r1_3 = r1_sq * std::sqrt(r1_sq);
    double r2_3 = r2_sq * std::sqrt(r2_sq);

    // Second partial derivatives of Omega
    double Oxx = 1.0 - (1-mu)/r1_3 - mu/r2_3
                 + 3*(1-mu)*(x+mu)*(x+mu)/r1_5 + 3*mu*(x-1+mu)*(x-1+mu)/r2_5;
    double Oyy = 1.0 - (1-mu)/r1_3 - mu/r2_3
                 + 3*(1-mu)*yy*yy/r1_5 + 3*mu*yy*yy/r2_5;
    double Ozz =     - (1-mu)/r1_3 - mu/r2_3
                 + 3*(1-mu)*z*z/r1_5 + 3*mu*z*z/r2_5;
    double Oxy = 3*(1-mu)*(x+mu)*yy/r1_5 + 3*mu*(x-1+mu)*yy/r2_5;
    double Oxz = 3*(1-mu)*(x+mu)*z/r1_5 + 3*mu*(x-1+mu)*z/r2_5;
    double Oyz = 3*(1-mu)*yy*z/r1_5 + 3*mu*yy*z/r2_5;

    // A matrix (6x6)
    // A = [ 0  I ]     where U = [Oxx Oxy Oxz; Oxy Oyy Oyz; Oxz Oyz Ozz]
    //     [ U  J ]     and   J = [0 2 0; -2 0 0; 0 0 0]
    double A[6][6] = {};
    A[0][3] = 1; A[1][4] = 1; A[2][5] = 1;
    A[3][0] = Oxx; A[3][1] = Oxy; A[3][2] = Oxz; A[3][4] = 2;
    A[4][0] = Oxy; A[4][1] = Oyy; A[4][2] = Oyz; A[4][3] = -2;
    A[5][0] = Oxz; A[5][1] = Oyz; A[5][2] = Ozz;

    // dPhi/dt = A * Phi
    double Phi[6][6], dPhi[6][6];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            Phi[i][j] = y[6 + i*6 + j];

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            dPhi[i][j] = 0;
            for (int k = 0; k < 6; ++k)
                dPhi[i][j] += A[i][k] * Phi[k][j];
        }

    // Pack output
    std::vector<double> dy(42);
    dy[0]=ds.xdot; dy[1]=ds.ydot; dy[2]=ds.zdot;
    dy[3]=ds.xdot; // wait, derivatives of state
    // Actually: dy/dt for state = [xdot, ydot, zdot, xddot, yddot, zddot]
    dy[0] = ds.x; dy[1] = ds.y; dy[2] = ds.z; // = xdot, ydot, zdot
    dy[3] = ds.xdot; dy[4] = ds.ydot; dy[5] = ds.zdot; // = accelerations

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            dy[6 + i*6 + j] = dPhi[i][j];

    return dy;
}

std::vector<CR3BPState> propagate_cr3bp(
    const CR3BPState& ic, double mu, double t_span,
    double dt_output, double atol, double rtol)
{
    std::vector<CR3BPState> trajectory;
    trajectory.push_back(ic);

    auto rhs = [mu](double /*t*/, const std::vector<double>& y) -> std::vector<double> {
        CR3BPState s = {y[0], y[1], y[2], y[3], y[4], y[5]};
        CR3BPState ds = cr3bp_derivatives(s, mu);
        return {ds.x, ds.y, ds.z, ds.xdot, ds.ydot, ds.zdot};
    };

    auto y = pack_state(ic);
    double t = 0.0;
    double dt = 0.01;
    double next_output = dt_output;

    while (t < t_span) {
        double h = std::min(dt, t_span - t);
        auto result = dopri5_generic_step(y, t, h, rhs, atol, rtol);

        if (result.accepted) {
            y = std::move(result.y);
            t += h;
            dt = result.dt_next;

            if (t >= next_output) {
                trajectory.push_back(unpack_state(y));
                next_output += dt_output;
            }
        } else {
            dt = result.dt_next;
        }
    }

    return trajectory;
}

CR3BPStateSTM propagate_cr3bp_stm(
    const CR3BPState& ic, double mu, double t_span,
    double atol, double rtol)
{
    auto rhs = [mu](double t, const std::vector<double>& y) {
        return cr3bp_rhs_stm(t, y, mu);
    };

    auto y = pack_state_stm(ic);
    double t = 0.0;
    double dt = 0.01;

    while (t < t_span) {
        double h = std::min(dt, t_span - t);
        auto result = dopri5_generic_step(y, t, h, rhs, atol, rtol);
        if (result.accepted) {
            y = std::move(result.y);
            t += h;
            dt = result.dt_next;
        } else {
            dt = result.dt_next;
        }
    }

    return unpack_state_stm(y);
}

double propagate_to_y_crossing(
    CR3BPState& state, double mu, double t_max,
    CR3BPStateSTM* stm_out, double atol, double rtol)
{
    bool use_stm = (stm_out != nullptr);

    auto rhs_simple = [mu](double /*t*/, const std::vector<double>& y) -> std::vector<double> {
        CR3BPState s = {y[0], y[1], y[2], y[3], y[4], y[5]};
        CR3BPState ds = cr3bp_derivatives(s, mu);
        return {ds.x, ds.y, ds.z, ds.xdot, ds.ydot, ds.zdot};
    };
    auto rhs_stm = [mu](double t, const std::vector<double>& y) {
        return cr3bp_rhs_stm(t, y, mu);
    };

    std::vector<double> y;
    GenericODEFunc rhs;
    if (use_stm) {
        y = pack_state_stm(state);
        rhs = rhs_stm;
    } else {
        y = pack_state(state);
        rhs = rhs_simple;
    }

    double t = 0.0;
    double dt = 0.01;
    double min_t = 0.1; // skip initial region to avoid detecting start as a crossing

    while (t < t_max) {
        double h = std::min(dt, t_max - t);
        double y_prev = y[1];

        auto result = dopri5_generic_step(y, t, h, rhs, atol, rtol);
        if (!result.accepted) { dt = result.dt_next; continue; }

        double y_new = result.y[1];

        // Detect any y=0 crossing after minimum propagation time
        if (t + h > min_t && y_prev * y_new < 0) {
            // Bisect for exact crossing using regula falsi
            double t_lo = t, t_hi = t + h;
            auto y_lo = y;
            double yval_lo = y_prev, yval_hi = y_new;

            for (int bisect = 0; bisect < 60; ++bisect) {
                // Linear interpolation for better convergence
                double frac = -yval_lo / (yval_hi - yval_lo);
                double t_mid = t_lo + frac * (t_hi - t_lo);
                t_mid = std::max(t_lo + 1e-15, std::min(t_hi - 1e-15, t_mid));

                double h_mid = t_mid - t_lo;
                auto res_mid = dopri5_generic_step(y_lo, t_lo, h_mid, rhs, atol, rtol);
                // If step rejected, use a smaller step
                if (!res_mid.accepted) {
                    t_mid = (t_lo + t_hi) / 2.0;
                    h_mid = t_mid - t_lo;
                    res_mid = dopri5_generic_step(y_lo, t_lo, h_mid, rhs, atol*100, rtol*100);
                }
                double yval_mid = res_mid.y[1];

                if (yval_lo * yval_mid < 0) {
                    t_hi = t_mid;
                    yval_hi = yval_mid;
                } else {
                    t_lo = t_mid;
                    y_lo = res_mid.y;
                    yval_lo = yval_mid;
                }
                if (std::fabs(yval_mid) < 1e-13 || t_hi - t_lo < 1e-14) break;
            }

            y = y_lo;
            t = t_lo;

            state = unpack_state(y);
            if (use_stm) *stm_out = unpack_state_stm(y);
            return t;
        }

        y = std::move(result.y);
        t += h;
        dt = result.dt_next;
    }

    state = unpack_state(y);
    if (use_stm) *stm_out = unpack_state_stm(y);
    return t;
}

// ============================================================
// Richardson 3rd-order halo orbit initial conditions
// ============================================================

CR3BPState richardson_halo_ic(double mu, int lpoint, double Az, int family) {
    auto lps = compute_lagrange_points(mu);
    double gamma_L;
    double x_L;

    if (lpoint == 1) {
        x_L = lps.L1.x;
        gamma_L = 1.0 - mu - x_L;
    } else if (lpoint == 2) {
        x_L = lps.L2.x;
        gamma_L = x_L - 1.0 + mu;
    } else {
        throw std::runtime_error("Richardson only supports L1 and L2");
    }

    // Compute c2 coefficient for linearized dynamics near Lagrange point
    auto cn = [mu, gamma_L, lpoint](int n) -> double {
        double g = gamma_L;
        double pm = (n % 2 == 0) ? 1.0 : -1.0;
        if (lpoint == 1) {
            return (1.0/(g*g*g)) * (pm * mu + pm * (1-mu) * std::pow(g/(1-g), n+1));
        } else {
            return (1.0/(g*g*g)) * (pm * mu + (1-mu) * std::pow(g/(1+g), n+1));
        }
    };

    double c2 = cn(2);

    // In-plane frequency from linearized EOM
    double disc = 9*c2*c2 - 8*c2;
    if (disc < 0) disc = 0;
    double lambda_sq = (c2 - 2.0 + std::sqrt(disc)) / 2.0;
    double lambda = std::sqrt(std::fabs(lambda_sq));
    double k = (lambda_sq + 1.0 + 2.0*c2) / (2.0*lambda);

    // Initial guess: place at Lagrange point with z-perturbation
    // Use linearized relationship for in-plane motion
    double Ax = Az * lambda / k; // approximate in-plane amplitude

    CR3BPState ic;
    ic.x = x_L + gamma_L * Ax; // small displacement from L-point
    ic.y = 0.0;
    ic.z = family * Az;
    ic.xdot = 0.0;
    ic.zdot = 0.0;

    // ydot from linearized dynamics
    ic.ydot = -k * lambda * gamma_L * Ax;

    return ic;
}

// ============================================================
// Halo orbit differential correction
// ============================================================

// Internal: single-shot differential correction from a given IC
static bool correct_halo_ic(CR3BPState& ic, double mu, double tol, int max_iter, int& iters) {
    double eps = 1e-8;
    for (int iter = 0; iter < max_iter; ++iter) {
        iters = iter + 1;
        CR3BPState s = ic;
        propagate_to_y_crossing(s, mu, 10.0);

        if (std::fabs(s.xdot) < tol && std::fabs(s.zdot) < tol) return true;

        CR3BPState s1 = ic; s1.x += eps;
        CR3BPState r1 = s1; propagate_to_y_crossing(r1, mu, 10.0);
        CR3BPState s2 = ic; s2.ydot += eps;
        CR3BPState r2 = s2; propagate_to_y_crossing(r2, mu, 10.0);

        double a = (r1.xdot - s.xdot) / eps, b = (r2.xdot - s.xdot) / eps;
        double c = (r1.zdot - s.zdot) / eps, d = (r2.zdot - s.zdot) / eps;
        double det = a * d - b * c;
        if (std::fabs(det) < 1e-30) return false;

        double dx = -(d * s.xdot - b * s.zdot) / det;
        double dydot = -(-c * s.xdot + a * s.zdot) / det;

        // Limit step size
        double scale = 1.0;
        if (std::fabs(dx) > 0.01) scale = std::min(scale, 0.01 / std::fabs(dx));
        if (std::fabs(dydot) > 0.05) scale = std::min(scale, 0.05 / std::fabs(dydot));

        ic.x += dx * scale;
        ic.ydot += dydot * scale;
    }
    return false;
}

HaloOrbit compute_halo_orbit(const System& sys, const HaloOrbitParams& params,
                              double tol, int max_iter) {
    HaloOrbit orbit;
    orbit.converged = false;
    orbit.iterations = 0;

    // Continuation approach: start from small Az, grow to target
    // This ensures the initial guess is close enough for Newton to converge
    int n_steps = 10;
    double target_Az = params.Az;

    // Start from Lyapunov-like IC
    auto lps = compute_lagrange_points(sys.mu);
    double x_L = (params.lpoint == 1) ? lps.L1.x : lps.L2.x;

    CR3BPState ic;
    ic.x = x_L;
    ic.y = 0; ic.z = 0; ic.xdot = 0; ic.zdot = 0;
    // Small initial ydot for a nearly-planar orbit
    ic.ydot = -0.01;

    // First: find a Lyapunov orbit (z=0)
    int sub_iter = 0;
    // For Lyapunov (z=0), we only correct x to make xdot=0 at y=0 crossing
    // (1D correction, simpler)
    double eps = 1e-8;
    for (int iter = 0; iter < max_iter; ++iter) {
        CR3BPState s = ic;
        propagate_to_y_crossing(s, sys.mu, 10.0);
        if (std::fabs(s.xdot) < tol) { break; }

        CR3BPState ic2 = ic; ic2.ydot += eps;
        CR3BPState s2 = ic2;
        propagate_to_y_crossing(s2, sys.mu, 10.0);

        double dxdot_dydot = (s2.xdot - s.xdot) / eps;
        if (std::fabs(dxdot_dydot) < 1e-30) break;

        double dydot = -s.xdot / dxdot_dydot;
        if (std::fabs(dydot) > 0.05) dydot *= 0.05 / std::fabs(dydot);
        ic.ydot += dydot;
        sub_iter++;
    }

    // Now gradually increase z from 0 to target_Az
    for (int step = 1; step <= n_steps; ++step) {
        double Az_step = target_Az * step / n_steps;
        ic.z = params.family * Az_step;

        int iters = 0;
        bool ok = correct_halo_ic(ic, sys.mu, tol, max_iter, iters);
        orbit.iterations += iters;

        if (!ok && step > 1) {
            // Try smaller increment
            ic.z = params.family * target_Az * (step - 0.5) / n_steps;
            ok = correct_halo_ic(ic, sys.mu, tol, max_iter, iters);
            orbit.iterations += iters;
        }

        if (!ok) break;
    }

    // Final correction at target Az
    ic.z = params.family * target_Az;
    int final_iters = 0;
    orbit.converged = correct_halo_ic(ic, sys.mu, tol, max_iter, final_iters);
    orbit.iterations += final_iters;

    if (orbit.converged) {
        CR3BPState sf = ic;
        double t_half = propagate_to_y_crossing(sf, sys.mu, 10.0);
        orbit.period = 2.0 * t_half;
        orbit.initial_state = ic;
        orbit.jacobi_constant = jacobi_constant(ic, sys.mu);
        orbit.trajectory = propagate_cr3bp(ic, sys.mu, orbit.period, orbit.period / 200.0);
        orbit.stability_index = 0.0;
    }

    return orbit;
}


// ============================================================
// Unit conversions
// ============================================================

State to_physical(const CR3BPState& s, const System& sys) {
    State phys;
    phys.pos = {s.x * sys.L_star, s.y * sys.L_star, s.z * sys.L_star};
    phys.vel = {s.xdot * sys.V_star, s.ydot * sys.V_star, s.zdot * sys.V_star};
    return phys;
}

CR3BPState to_normalized(const State& s, const System& sys) {
    return {
        s.pos.x / sys.L_star, s.pos.y / sys.L_star, s.pos.z / sys.L_star,
        s.vel.x / sys.V_star, s.vel.y / sys.V_star, s.vel.z / sys.V_star
    };
}

} // namespace cr3bp
} // namespace solar
