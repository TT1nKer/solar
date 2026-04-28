#include "solar/kepler.h"
#include "solar/constants.h"
#include <cmath>
#include <stdexcept>

namespace solar {

double solve_kepler(double M, double e, double tol, int max_iter) {
    // Normalize M to [0, 2*pi)
    M = std::fmod(M, constants::TWO_PI);
    if (M < 0.0) M += constants::TWO_PI;

    // Initial guess
    double E = (e > 0.8) ? constants::PI : M;

    // Newton-Raphson iteration: f(E) = E - e*sin(E) - M
    for (int i = 0; i < max_iter; ++i) {
        double f  = E - e * std::sin(E) - M;
        double fp = 1.0 - e * std::cos(E);
        double dE = f / fp;
        E -= dE;
        if (std::fabs(dE) < tol) return E;
    }

    return E; // best estimate if not converged
}

double mean_to_true_anomaly(double M, double e) {
    double E = solve_kepler(M, e);

    // True anomaly from eccentric anomaly
    double half_E = E / 2.0;
    double nu = 2.0 * std::atan2(
        std::sqrt(1.0 + e) * std::sin(half_E),
        std::sqrt(1.0 - e) * std::cos(half_E)
    );

    return nu;
}

double true_to_eccentric_anomaly(double nu, double e) {
    return 2.0 * std::atan2(
        std::sqrt(1.0 - e) * std::sin(nu / 2.0),
        std::sqrt(1.0 + e) * std::cos(nu / 2.0)
    );
}

double eccentric_to_mean_anomaly(double E, double e) {
    return E - e * std::sin(E);
}

State elements_to_state(const OrbitalElements& elem, double mu, double true_anomaly) {
    double a = elem.a;
    double e = elem.e;
    double i = elem.i;
    double omega = elem.omega;
    double Omega = elem.Omega;

    // Distance from focus
    double p = a * (1.0 - e * e); // semi-latus rectum
    double r = p / (1.0 + e * std::cos(true_anomaly));

    // Position in perifocal frame (PQW)
    double cos_nu = std::cos(true_anomaly);
    double sin_nu = std::sin(true_anomaly);

    double x_pf = r * cos_nu;
    double y_pf = r * sin_nu;

    // Velocity in perifocal frame
    double sqrt_mu_p = std::sqrt(mu / p);
    double vx_pf = -sqrt_mu_p * sin_nu;
    double vy_pf =  sqrt_mu_p * (e + cos_nu);

    // Rotation from perifocal to heliocentric ecliptic frame
    double cos_O = std::cos(Omega);
    double sin_O = std::sin(Omega);
    double cos_w = std::cos(omega);
    double sin_w = std::sin(omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);

    // Rotation matrix columns (perifocal P and Q directions in ecliptic frame)
    double Px = cos_O * cos_w - sin_O * sin_w * cos_i;
    double Py = sin_O * cos_w + cos_O * sin_w * cos_i;
    double Pz = sin_w * sin_i;

    double Qx = -cos_O * sin_w - sin_O * cos_w * cos_i;
    double Qy = -sin_O * sin_w + cos_O * cos_w * cos_i;
    double Qz = cos_w * sin_i;

    State s;
    s.pos.x = Px * x_pf + Qx * y_pf;
    s.pos.y = Py * x_pf + Qy * y_pf;
    s.pos.z = Pz * x_pf + Qz * y_pf;

    s.vel.x = Px * vx_pf + Qx * vy_pf;
    s.vel.y = Py * vx_pf + Qy * vy_pf;
    s.vel.z = Pz * vx_pf + Qz * vy_pf;

    return s;
}

State elements_to_state_at(const OrbitalElements& elem, double mu, double jd) {
    double n = mean_motion(elem.a, mu);
    double dt = (jd - elem.epoch) * constants::DAY; // seconds since epoch
    double M = elem.M0 + n * dt;
    double nu = mean_to_true_anomaly(M, elem.e);
    return elements_to_state(elem, mu, nu);
}

OrbitalElements state_to_elements(const State& s, double mu, double epoch) {
    Vec3 r = s.pos;
    Vec3 v = s.vel;

    double r_mag = r.norm();
    double v_mag = v.norm();

    // Specific angular momentum
    Vec3 h = r.cross(v);
    double h_mag = h.norm();

    // Node vector (z-hat cross h)
    Vec3 n = {-h.y, h.x, 0.0};
    double n_mag = n.norm();

    // Eccentricity vector
    Vec3 e_vec = (v.cross(h) / mu) - r.normalized();
    double e = e_vec.norm();

    // Semi-major axis (vis-viva)
    double energy = 0.5 * v_mag * v_mag - mu / r_mag;
    double a = -mu / (2.0 * energy);

    // Inclination
    double i = std::acos(h.z / h_mag);

    // Longitude of ascending node
    double Omega = 0.0;
    if (n_mag > 1e-10) {
        Omega = std::acos(n.x / n_mag);
        if (n.y < 0.0) Omega = constants::TWO_PI - Omega;
    }

    // Argument of periapsis
    double omega = 0.0;
    if (n_mag > 1e-10 && e > 1e-10) {
        omega = std::acos(n.dot(e_vec) / (n_mag * e));
        if (e_vec.z < 0.0) omega = constants::TWO_PI - omega;
    }

    // True anomaly
    double nu = 0.0;
    if (e > 1e-10) {
        nu = std::acos(e_vec.dot(r) / (e * r_mag));
        if (r.dot(v) < 0.0) nu = constants::TWO_PI - nu;
    }

    // Mean anomaly
    double E = true_to_eccentric_anomaly(nu, e);
    double M0 = eccentric_to_mean_anomaly(E, e);

    return {a, e, i, omega, Omega, M0, epoch};
}

double orbital_period(double a, double mu) {
    return constants::TWO_PI * std::sqrt(a * a * a / mu);
}

double mean_motion(double a, double mu) {
    return std::sqrt(mu / (a * a * a));
}

double vis_viva_speed(double r, double a, double mu) {
    return std::sqrt(mu * (2.0 / r - 1.0 / a));
}

} // namespace solar
