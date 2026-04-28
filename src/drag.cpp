#include "solar/drag.h"
#include "solar/atmosphere.h"
#include "solar/iau_rotation.h"
#include "solar/jpl_de.h"
#include "solar/constants.h"
#include <cmath>

namespace solar {

AtmosphericDrag::AtmosphericDrag(std::vector<DragTarget> targets)
    : targets_(std::move(targets)) {}

void AtmosphericDrag::compute(
    const std::vector<Body>& bodies,
    double time,
    std::vector<Vec3>& acc) const
{
    for (const auto& target : targets_) {
        int si = target.body_index;
        if (si < 0 || si >= static_cast<int>(bodies.size())) continue;

        // Find parent body (the one with the atmosphere)
        int pi = bodies[si].parent_index;
        if (pi < 0 || pi >= static_cast<int>(bodies.size())) continue;

        // Get atmosphere profile
        const AtmosphereProfile* atm = atmosphere_for_body(bodies[pi].name);
        if (!atm) continue;

        // Position relative to parent
        Vec3 r_rel = bodies[si].state.pos - bodies[pi].state.pos;
        double r_mag = r_rel.norm();
        double altitude = r_mag - atm->surface_radius;

        // Compute density
        double rho = atmospheric_density(*atm, altitude);
        if (rho <= 0.0) continue;

        // Velocity relative to parent
        Vec3 v_rel = bodies[si].state.vel - bodies[pi].state.vel;

        // Account for atmosphere co-rotation
        // omega is the angular velocity of the parent body
        int iau_id = body_name_to_de_id(bodies[pi].name);
        if (iau_id > 0) {
            Vec3 omega;
            iau_body_fixed_matrix(iau_id, time / constants::DAY + constants::J2000, omega);
            // Atmosphere rotates with the body: v_atm = omega x r_rel
            Vec3 v_atm = omega.cross(r_rel);
            v_rel = v_rel - v_atm;
        }

        double v_rel_mag = v_rel.norm();
        if (v_rel_mag < 1e-10) continue;

        // Drag acceleration: a = -0.5 * rho * v^2 * Cd * A / m * v_hat
        // Units: rho (kg/m^3), v (km/s), A (m^2), m (kg)
        // a = -0.5 * rho * (v*1000)^2 * Cd * A / m / 1000 [km/s^2]
        // Simplifies to: a = -0.5 * rho * v^2 * Cd * A / m * 1000
        double mass_kg = bodies[si].mass;
        if (mass_kg <= 0.0) continue;

        double a_mag = 0.5 * rho * v_rel_mag * v_rel_mag * target.Cd
                     * target.area_m2 / mass_kg * 1000.0;

        Vec3 v_hat = v_rel / v_rel_mag;
        acc[si] -= v_hat * a_mag;
    }
}

} // namespace solar
