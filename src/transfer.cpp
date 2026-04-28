#include "solar/transfer.h"
#include "solar/kepler.h"
#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <cmath>
#include <stdexcept>

namespace solar {

TransferResult hohmann_transfer(double r1, double r2, double mu) {
    TransferResult result;

    // Transfer orbit semi-major axis
    result.a_transfer = (r1 + r2) / 2.0;

    // Eccentricity of transfer orbit
    result.e_transfer = std::fabs(r2 - r1) / (r1 + r2);

    // Velocities on circular orbits
    double v_circ1 = std::sqrt(mu / r1);
    double v_circ2 = std::sqrt(mu / r2);

    // Velocities on transfer orbit at departure and arrival
    double v_trans1 = vis_viva_speed(r1, result.a_transfer, mu);
    double v_trans2 = vis_viva_speed(r2, result.a_transfer, mu);

    // Delta-v at each burn
    result.dv1 = std::fabs(v_trans1 - v_circ1);
    result.dv2 = std::fabs(v_circ2 - v_trans2);
    result.dv_total = result.dv1 + result.dv2;

    // Time of flight: half the transfer orbit period
    result.tof = constants::PI * std::sqrt(result.a_transfer * result.a_transfer * result.a_transfer / mu);

    return result;
}

TransferResult bielliptic_transfer(double r1, double r2, double r_intermediate, double mu) {
    TransferResult result;

    // First transfer ellipse: r1 to r_intermediate
    double a1 = (r1 + r_intermediate) / 2.0;
    // Second transfer ellipse: r_intermediate to r2
    double a2 = (r2 + r_intermediate) / 2.0;

    double v_circ1 = std::sqrt(mu / r1);
    double v_circ2 = std::sqrt(mu / r2);

    double v_trans1_dep = vis_viva_speed(r1, a1, mu);
    double v_trans1_arr = vis_viva_speed(r_intermediate, a1, mu);
    double v_trans2_dep = vis_viva_speed(r_intermediate, a2, mu);
    double v_trans2_arr = vis_viva_speed(r2, a2, mu);

    result.dv1 = std::fabs(v_trans1_dep - v_circ1);
    double dv_intermediate = std::fabs(v_trans2_dep - v_trans1_arr);
    result.dv2 = std::fabs(v_circ2 - v_trans2_arr);
    result.dv_total = result.dv1 + dv_intermediate + result.dv2;

    // Total TOF: half-periods of both transfer ellipses
    result.tof = constants::PI * std::sqrt(a1 * a1 * a1 / mu)
               + constants::PI * std::sqrt(a2 * a2 * a2 / mu);

    result.a_transfer = a1; // first ellipse semi-major axis
    result.e_transfer = (r_intermediate - r1) / (r_intermediate + r1);

    return result;
}

// Stumpff functions for universal variable formulation
static double stumpff_c2(double psi) {
    if (std::fabs(psi) < 1e-6) {
        return 1.0 / 2.0 - psi / 24.0 + psi * psi / 720.0;
    }
    if (psi > 0) {
        double sp = std::sqrt(psi);
        return (1.0 - std::cos(sp)) / psi;
    } else {
        double sp = std::sqrt(-psi);
        return (1.0 - std::cosh(sp)) / psi;
    }
}

static double stumpff_c3(double psi) {
    if (std::fabs(psi) < 1e-6) {
        return 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
    }
    if (psi > 0) {
        double sp = std::sqrt(psi);
        return (sp - std::sin(sp)) / (psi * sp);
    } else {
        double sp = std::sqrt(-psi);
        return (std::sinh(sp) - sp) / (-psi * sp);
    }
}

LambertResult solve_lambert(const Vec3& r1, const Vec3& r2, double tof, double mu,
                            bool prograde) {
    // Universal variable Lambert solver (Curtis, "Orbital Mechanics for Engineering Students")
    LambertResult result;
    result.converged = false;
    result.iterations = 0;

    double r1_mag = r1.norm();
    double r2_mag = r2.norm();

    // Determine transfer angle
    double cos_dnu = r1.dot(r2) / (r1_mag * r2_mag);
    cos_dnu = std::max(-1.0, std::min(1.0, cos_dnu));

    // Cross product z-component determines prograde/retrograde
    double cross_z = r1.x * r2.y - r1.y * r2.x;

    double A;
    if (prograde) {
        A = (cross_z >= 0) ? 1.0 : -1.0;
    } else {
        A = (cross_z < 0) ? 1.0 : -1.0;
    }
    A *= std::sqrt(r1_mag * r2_mag * (1.0 + cos_dnu));

    if (std::fabs(A) < 1e-14) {
        return result; // degenerate (180-degree transfer)
    }

    // Helper: compute y(z) and F(z) = [y/c2]^1.5 * c3 + A*sqrt(y) - sqrt(mu)*tof
    double sqrt_mu = std::sqrt(mu);

    // Use bisection to find robust initial bounds, then refine
    // F(z) is monotonically decreasing in z for single-revolution transfers
    // Find z such that F(z) = 0

    // Start with bisection on a wide range
    double z_low = -4.0 * constants::PI * constants::PI;
    double z_up  = 4.0 * constants::PI * constants::PI;

    // Ensure y > 0 at lower bound
    {
        double c2 = stumpff_c2(z_low);
        double c3 = stumpff_c3(z_low);
        double y = r1_mag + r2_mag + A * (z_low * c3 - 1.0) / std::sqrt(c2);
        while (y < 0.0) {
            z_low *= 0.5;
            c2 = stumpff_c2(z_low);
            c3 = stumpff_c3(z_low);
            y = r1_mag + r2_mag + A * (z_low * c3 - 1.0) / std::sqrt(c2);
        }
    }

    int max_iter = 200;
    double tol = 1e-8;

    // Evaluate F at bounds to ensure they bracket the root
    auto eval_F = [&](double z) -> double {
        double c2 = stumpff_c2(z);
        double c3 = stumpff_c3(z);
        double y = r1_mag + r2_mag + A * (z * c3 - 1.0) / std::sqrt(c2);
        if (y < 0.0) return -1e30; // signal invalid
        double x = std::sqrt(y / c2);
        return x * x * x * c3 + A * std::sqrt(y) - sqrt_mu * tof;
    };

    double F_low = eval_F(z_low);
    double F_up = eval_F(z_up);

    // If F doesn't change sign, try wider bounds
    if (F_low * F_up > 0) {
        // Try extending
        for (int i = 0; i < 20; ++i) {
            z_up *= 2.0;
            F_up = eval_F(z_up);
            if (F_low * F_up <= 0) break;
        }
    }

    if (F_low * F_up > 0) {
        return result; // cannot bracket root
    }

    // Bisection method (very robust)
    double z = (z_low + z_up) / 2.0;
    for (int iter = 0; iter < max_iter; ++iter) {
        result.iterations = iter + 1;

        double c2 = stumpff_c2(z);
        double c3 = stumpff_c3(z);
        double y = r1_mag + r2_mag + A * (z * c3 - 1.0) / std::sqrt(c2);

        if (y < 0.0) {
            z_low = z;
            z = (z_low + z_up) / 2.0;
            continue;
        }

        double x = std::sqrt(y / c2);
        double F = x * x * x * c3 + A * std::sqrt(y) - sqrt_mu * tof;

        if (std::fabs(F) < tol) {
            result.converged = true;
            break;
        }

        if (F < 0) {
            z_low = z;
        } else {
            z_up = z;
        }

        z = (z_low + z_up) / 2.0;

        if (std::fabs(z_up - z_low) < 1e-14) {
            result.converged = true;
            break;
        }
    }

    // Compute velocities from converged z
    double c2 = stumpff_c2(z);
    double c3 = stumpff_c3(z);
    double y = r1_mag + r2_mag + A * (z * c3 - 1.0) / std::sqrt(c2);

    double f = 1.0 - y / r1_mag;
    double g = A * std::sqrt(y / mu);
    double g_dot = 1.0 - y / r2_mag;

    result.v1 = (r2 - r1 * f) / g;
    result.v2 = (r2 * g_dot - r1) / g;

    return result;
}

TransferResult plan_transfer(const Body& departure, const Body& arrival,
                             double jd_depart, double jd_arrive, double mu_central) {
    // Get positions at departure and arrival dates
    State s1 = compute_position(departure, jd_depart);
    State s2 = compute_position(arrival, jd_arrive);

    double tof = (jd_arrive - jd_depart) * constants::DAY; // seconds

    // Solve Lambert problem
    LambertResult lambert = solve_lambert(s1.pos, s2.pos, tof, mu_central);

    TransferResult result;
    if (!lambert.converged) {
        result.dv1 = -1.0;
        result.dv2 = -1.0;
        result.dv_total = -1.0;
        result.tof = tof;
        result.a_transfer = 0.0;
        result.e_transfer = 0.0;
        return result;
    }

    // Delta-v at departure: difference between transfer orbit and planet's velocity
    Vec3 dv1_vec = lambert.v1 - s1.vel;
    Vec3 dv2_vec = s2.vel - lambert.v2;

    result.dv1 = dv1_vec.norm();
    result.dv2 = dv2_vec.norm();
    result.dv_total = result.dv1 + result.dv2;
    result.tof = tof;

    // Compute transfer orbit elements
    State transfer_state = {s1.pos, lambert.v1};
    OrbitalElements transfer_elem = state_to_elements(transfer_state, mu_central);
    result.a_transfer = transfer_elem.a;
    result.e_transfer = transfer_elem.e;

    return result;
}

} // namespace solar
