#include "solar/j2.h"
#include "solar/ephemeris.h"
#include <cmath>

namespace solar {

J2Perturbation::J2Perturbation(std::vector<J2Entry> entries)
    : entries_(std::move(entries)) {}

J2Perturbation J2Perturbation::solar_system_defaults(const std::vector<Body>& bodies) {
    std::vector<J2Entry> entries;

    auto add = [&](const char* name, double j2, double r_eq) {
        const Body* b = find_body(bodies, name);
        if (b) {
            int idx = static_cast<int>(b - bodies.data());
            entries.push_back({idx, j2, r_eq});
        }
    };

    // Major solar system J2 values
    add("Sun",     2.0e-7,     695700.0);
    add("Earth",   1.08263e-3, 6378.137);
    add("Mars",    1.96045e-3, 3396.2);
    add("Jupiter", 1.4736e-2,  71492.0);
    add("Saturn",  1.6298e-2,  60268.0);
    add("Uranus",  3.343e-3,   25559.0);
    add("Neptune", 3.411e-3,   24764.0);

    return J2Perturbation(std::move(entries));
}

void J2Perturbation::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const
{
    size_t n = bodies.size();

    for (const auto& entry : entries_) {
        int j = entry.body_index;
        if (j < 0 || j >= static_cast<int>(n)) continue;

        double mu = bodies[j].mu;
        double J2 = entry.J2;
        double Re = entry.R_eq;
        double Re2 = Re * Re;

        // Apply J2 perturbation from body j to all other bodies
        for (size_t i = 0; i < n; ++i) {
            if (static_cast<int>(i) == j) continue;

            Vec3 r = bodies[i].state.pos - bodies[j].state.pos;
            double r_mag = r.norm();
            double r2 = r_mag * r_mag;
            double r5 = r2 * r2 * r_mag;

            // z-component relative to oblate body (assume pole aligned with ecliptic z)
            double z = r.z;
            double z2_r2 = (z * z) / r2;

            // J2 acceleration components
            // a = -3/2 * J2 * mu * Re^2 / r^5 * [(1 - 5*z^2/r^2)*r + 2*z*z_hat]
            double coeff = -1.5 * J2 * mu * Re2 / r5;

            double factor_xy = 1.0 - 5.0 * z2_r2;
            acc[i].x += coeff * factor_xy * r.x;
            acc[i].y += coeff * factor_xy * r.y;
            acc[i].z += coeff * (3.0 - 5.0 * z2_r2) * r.z;
        }
    }
}

} // namespace solar
