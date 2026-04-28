#include "solar/srp.h"
#include "solar/constants.h"
#include <cmath>

namespace solar {

// Solar constant at 1 AU (W/m^2)
static constexpr double P_SUN_1AU = 1361.0;

SolarRadiationPressure::SolarRadiationPressure(
    int sun_index, std::vector<SRPTarget> targets, bool enable_shadow)
    : sun_index_(sun_index), targets_(std::move(targets)),
      enable_shadow_(enable_shadow) {}

double SolarRadiationPressure::shadow_factor(
    const Vec3& sc_pos, const Vec3& sun_pos,
    const std::vector<Body>& bodies) const
{
    // Cylindrical shadow model: check if any body occludes the Sun
    Vec3 sc_to_sun = sun_pos - sc_pos;
    double dist_to_sun = sc_to_sun.norm();
    Vec3 sun_dir = sc_to_sun / dist_to_sun;

    for (size_t i = 0; i < bodies.size(); ++i) {
        if (static_cast<int>(i) == sun_index_) continue;
        if (bodies[i].type == BodyType::Spacecraft) continue;
        if (bodies[i].properties.radius_km <= 0.0) continue;

        Vec3 sc_to_body = bodies[i].state.pos - sc_pos;
        double proj = sc_to_body.dot(sun_dir);

        // Body must be between spacecraft and Sun
        if (proj <= 0.0 || proj >= dist_to_sun) continue;

        // Perpendicular distance from Sun-line to body center
        Vec3 closest = sc_pos + sun_dir * proj;
        double perp_dist = (bodies[i].state.pos - closest).norm();

        if (perp_dist < bodies[i].properties.radius_km) {
            return 0.0; // in shadow
        }
    }

    return 1.0; // in sunlight
}

void SolarRadiationPressure::compute(
    const std::vector<Body>& bodies,
    double /*time*/,
    std::vector<Vec3>& acc) const
{
    if (sun_index_ < 0 || sun_index_ >= static_cast<int>(bodies.size())) return;

    const Vec3& sun_pos = bodies[sun_index_].state.pos;

    for (const auto& target : targets_) {
        int si = target.body_index;
        if (si < 0 || si >= static_cast<int>(bodies.size())) continue;

        Vec3 r_sun = bodies[si].state.pos - sun_pos; // spacecraft relative to Sun
        double r_mag = r_sun.norm();
        if (r_mag < 1.0) continue; // avoid division by zero

        // Shadow check
        double shadow = enable_shadow_
            ? shadow_factor(bodies[si].state.pos, sun_pos, bodies) : 1.0;
        if (shadow == 0.0) continue;

        // SRP acceleration:
        // F/m = (P_sun / c) * Cr * A / m * (AU/r)^2, directed away from Sun
        // P_sun at distance r: P = P_SUN_1AU * (AU/r)^2
        // Radiation pressure: P/c where c in m/s
        // Convert to km/s^2: divide by 1000

        double AU_m = constants::AU * 1000.0; // AU in meters
        double r_m = r_mag * 1000.0;          // distance in meters
        double c_m = constants::C_LIGHT * 1000.0; // speed of light in m/s

        double P_at_r = P_SUN_1AU * (AU_m * AU_m) / (r_m * r_m); // W/m^2 at distance
        double pressure = P_at_r / c_m; // N/m^2

        double mass_kg = bodies[si].mass;
        if (mass_kg <= 0.0) continue;

        // Acceleration in m/s^2, then convert to km/s^2
        double a_mag = shadow * pressure * target.Cr * target.area_m2 / mass_kg / 1000.0;

        Vec3 r_hat = r_sun / r_mag;
        acc[si] += r_hat * a_mag; // directed away from Sun
    }
}

} // namespace solar
