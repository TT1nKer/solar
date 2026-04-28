#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace solar;
using namespace solar::constants;

int main() {
    auto bodies = load_solar_system();

    // J2000.0 = 2000-01-01 12:00 TT = JD 2451545.0
    double jd = J2000;

    std::cout << std::scientific << std::setprecision(10);
    std::cout << "=== Our simulator at J2000.0 (JD " << std::fixed << std::setprecision(1) << jd << ") ===\n\n";

    const char* names[] = {"Earth", "Mars"};
    for (const char* name : names) {
        const Body* b = find_body(bodies, name);
        if (!b) continue;
        State s = compute_position(*b, jd);
        std::cout << name << ":\n";
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  X = " << s.pos.x << " km\n";
        std::cout << "  Y = " << s.pos.y << " km\n";
        std::cout << "  Z = " << s.pos.z << " km\n";
        std::cout << "  VX = " << s.vel.x << " km/s\n";
        std::cout << "  VY = " << s.vel.y << " km/s\n";
        std::cout << "  VZ = " << s.vel.z << " km/s\n";
        std::cout << "  |r| = " << s.pos.norm() << " km (" << s.pos.norm() / AU << " AU)\n";
        std::cout << "  |v| = " << s.vel.norm() << " km/s\n\n";
    }

    // JPL Horizons values (ground truth)
    std::cout << "=== JPL Horizons (ground truth) ===\n\n";

    std::cout << "Earth:\n";
    std::cout << "  X = -2.6499034e+07 km\n";
    std::cout << "  Y =  1.4469730e+08 km\n";
    std::cout << "  Z = -6.1115e+02 km\n";
    std::cout << "  VX = -2.9794260e+01 km/s\n";
    std::cout << "  VY = -5.4692949e+00 km/s\n";
    std::cout << "  VZ =  1.8178e-04 km/s\n\n";

    std::cout << "Mars:\n";
    std::cout << "  X =  2.0804814e+08 km\n";
    std::cout << "  Y = -2.0070526e+06 km\n";
    std::cout << "  Z = -5.1562890e+06 km\n";
    std::cout << "  VX =  1.1626724e+00 km/s\n";
    std::cout << "  VY =  2.6296065e+01 km/s\n";
    std::cout << "  VZ =  5.2229702e-01 km/s\n\n";

    // Compute errors
    std::cout << "=== Comparison ===\n\n";

    struct Ref { double x, y, z, vx, vy, vz; };
    Ref earth_jpl = {-2.6499034e7, 1.4469730e8, -6.1115e2, -2.9794260e1, -5.4692949e0, 1.8178e-4};
    Ref mars_jpl  = {2.0804814e8, -2.0070526e6, -5.1562890e6, 1.1626724e0, 2.6296065e1, 5.2229702e-1};

    auto compare = [](const char* name, const State& ours, const Ref& jpl) {
        double pos_err = std::sqrt(
            (ours.pos.x - jpl.x) * (ours.pos.x - jpl.x) +
            (ours.pos.y - jpl.y) * (ours.pos.y - jpl.y) +
            (ours.pos.z - jpl.z) * (ours.pos.z - jpl.z));
        double r_jpl = std::sqrt(jpl.x*jpl.x + jpl.y*jpl.y + jpl.z*jpl.z);
        double vel_err = std::sqrt(
            (ours.vel.x - jpl.vx) * (ours.vel.x - jpl.vx) +
            (ours.vel.y - jpl.vy) * (ours.vel.y - jpl.vy) +
            (ours.vel.z - jpl.vz) * (ours.vel.z - jpl.vz));
        double v_jpl = std::sqrt(jpl.vx*jpl.vx + jpl.vy*jpl.vy + jpl.vz*jpl.vz);

        std::cout << std::fixed;
        std::cout << name << ":\n";
        std::cout << "  Position error: " << std::setprecision(0) << pos_err << " km ("
                  << std::setprecision(4) << (pos_err / r_jpl) * 100.0 << "%)\n";
        std::cout << "  Velocity error: " << std::setprecision(4) << vel_err << " km/s ("
                  << std::setprecision(4) << (vel_err / v_jpl) * 100.0 << "%)\n\n";
    };

    const Body* earth = find_body(bodies, "Earth");
    const Body* mars = find_body(bodies, "Mars");
    compare("Earth", compute_position(*earth, jd), earth_jpl);
    compare("Mars", compute_position(*mars, jd), mars_jpl);

    return 0;
}
