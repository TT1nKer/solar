#include "solar/spacecraft.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>

namespace solar {

double Spacecraft::burn(double duration_s, const Vec3& direction, std::vector<Body>& bodies) {
    if (fuel_mass <= 0.0 || body_index < 0) return 0.0;

    double ve = exhaust_velocity_kms();           // km/s
    double mdot = thrust_N / (ve * 1e3) / 1000.0; // mass flow rate (kg/s), thrust_N / (ve_m/s)
    mdot = thrust_N / (isp * 9.80665);            // correct: mdot = F / (Isp * g0)

    double fuel_needed = mdot * duration_s;
    double actual_fuel = std::min(fuel_needed, fuel_mass);

    double m0 = total_mass();
    fuel_mass -= actual_fuel;
    double m1 = total_mass();

    // Tsiolkovsky: actual delta-v achieved
    double dv = ve * std::log(m0 / m1);

    // Apply velocity change to the Body
    Vec3 dir = direction.normalized();
    bodies[body_index].state.vel += dir * dv;
    bodies[body_index].mass = total_mass();

    return dv;
}

Spacecraft create_spacecraft(std::vector<Body>& bodies,
                             const std::string& name,
                             const State& initial_state,
                             int parent_index,
                             double dry_mass, double fuel_mass,
                             double thrust_N, double isp) {
    Body b;
    b.name = name;
    b.mass = dry_mass + fuel_mass;
    b.mu = 0.0;  // spacecraft has negligible gravitational pull
    b.elements = {0, 0, 0, 0, 0, 0, 0}; // not meaningful for spacecraft
    b.state = initial_state;
    b.type = BodyType::Spacecraft;
    b.parent_index = parent_index;
    if (parent_index >= 0 && parent_index < static_cast<int>(bodies.size())) {
        b.mu_parent = bodies[parent_index].mu;
    }

    int idx = static_cast<int>(bodies.size());
    bodies.push_back(b);

    Spacecraft sc;
    sc.name = name;
    sc.dry_mass = dry_mass;
    sc.fuel_mass = fuel_mass;
    sc.thrust_N = thrust_N;
    sc.isp = isp;
    sc.body_index = idx;

    return sc;
}

} // namespace solar
