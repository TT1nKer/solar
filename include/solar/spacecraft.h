#pragma once

#include "solar/body.h"
#include <vector>
#include <string>
#include <cmath>

namespace solar {

struct Spacecraft {
    std::string name;
    double dry_mass;      // kg (structural mass without fuel)
    double fuel_mass;     // kg (current propellant mass)
    double thrust_N;      // engine thrust in Newtons
    double isp;           // specific impulse in seconds
    int body_index;       // index of this craft's Body in the bodies vector

    double total_mass() const { return dry_mass + fuel_mass; }

    // Exhaust velocity (km/s)
    double exhaust_velocity_kms() const {
        return isp * 9.80665 / 1000.0;
    }

    // Maximum delta-v remaining (Tsiolkovsky rocket equation) in km/s
    double delta_v_remaining() const {
        if (fuel_mass <= 0.0) return 0.0;
        return exhaust_velocity_kms() * std::log(total_mass() / dry_mass);
    }

    // Apply a burn: reduces fuel, returns actual delta-v achieved (km/s)
    // duration_s: burn duration in seconds
    // direction: unit vector of thrust direction
    // Updates fuel_mass and the Body's mass in the bodies vector
    double burn(double duration_s, const Vec3& direction, std::vector<Body>& bodies);
};

// Create a spacecraft and insert its Body into the bodies vector.
// Returns the Spacecraft struct with body_index set.
// initial_state: position/velocity (in parent's reference frame)
// parent_index: which body this spacecraft orbits (-1 for heliocentric)
Spacecraft create_spacecraft(std::vector<Body>& bodies,
                             const std::string& name,
                             const State& initial_state,
                             int parent_index,
                             double dry_mass, double fuel_mass,
                             double thrust_N, double isp);

} // namespace solar
