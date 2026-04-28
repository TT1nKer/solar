#include "solar/moons.h"
#include "solar/kepler.h"
#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <algorithm>
#include <cctype>

namespace solar {

// Find index of a body by name (case-insensitive)
static int find_index(const std::vector<Body>& bodies, const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    for (size_t i = 0; i < bodies.size(); ++i) {
        std::string bn = bodies[i].name;
        std::transform(bn.begin(), bn.end(), bn.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (bn == lower) return static_cast<int>(i);
    }
    return -1;
}

void load_moons(std::vector<Body>& bodies) {
    using namespace constants;

    // Helper to add a moon
    // a_km: semi-major axis in km (relative to parent)
    // e: eccentricity
    // i_deg: inclination to parent's equatorial plane (deg)
    // Omega_deg: longitude of ascending node (deg)
    // omega_deg: argument of periapsis (deg)
    // M0_deg: mean anomaly at J2000 (deg)
    auto make_moon = [&](const std::string& name, const std::string& parent_name,
                         double mass_kg, double mu_val, double mu_par,
                         double a_km, double e, double i_deg,
                         double Omega_deg, double omega_deg, double M0_deg,
                         double radius, double surf_grav, double esc_vel,
                         bool atmo = false, double atm_press = 0.0) {
        int pi = find_index(bodies, parent_name);
        if (pi < 0) return;

        Body b;
        b.name = name;
        b.mass = mass_kg;
        b.mu = mu_val;
        b.type = BodyType::Moon;
        b.parent_index = pi;
        b.mu_parent = mu_par;

        b.elements.a = a_km;
        b.elements.e = e;
        b.elements.i = i_deg * DEG2RAD;
        b.elements.omega = omega_deg * DEG2RAD;
        b.elements.Omega = Omega_deg * DEG2RAD;
        b.elements.M0 = M0_deg * DEG2RAD;
        b.elements.epoch = J2000;

        // Compute initial state relative to parent at J2000
        b.state = elements_to_state_at(b.elements, mu_par, J2000);

        b.properties.radius_km = radius;
        b.properties.surface_gravity = surf_grav;
        b.properties.escape_velocity = esc_vel;
        b.properties.has_atmosphere = atmo;
        b.properties.atm_surface_pressure = atm_press;

        bodies.push_back(b);
    };

    // ==============================
    // Earth's Moon (Luna)
    // ==============================
    // Elements: relative to Earth, J2000 (ecliptic-referenced approximation)
    make_moon("Moon", "Earth",
              7.342e22, MU_MOON, MU_EARTH,
              384400.0,     // a (km)
              0.0549,       // e
              5.145,        // i (deg) — to ecliptic
              125.08,       // Omega (deg)
              318.15,       // omega (deg)
              134.96,       // M0 (deg) at J2000
              1737.4,       // radius (km)
              1.62,         // surface gravity (m/s^2)
              2.38,         // escape velocity (km/s)
              false, 0.0);

    // ==============================
    // Mars moons
    // ==============================
    make_moon("Phobos", "Mars",
              1.0659e16, MU_PHOBOS, MU_MARS,
              9376.0,       // a (km)
              0.0151,       // e
              1.093,        // i (deg) — to Mars equator
              164.93,       // Omega (deg)
              150.06,       // omega (deg)
              92.47,        // M0 (deg)
              11.267,       // radius (km)
              0.0057,       // surface gravity (m/s^2)
              0.01139,      // escape velocity (km/s)
              false, 0.0);

    make_moon("Deimos", "Mars",
              1.4762e15, MU_DEIMOS, MU_MARS,
              23463.2,      // a (km)
              0.00033,      // e
              0.93,         // i (deg)
              339.60,       // Omega (deg)
              290.50,       // omega (deg)
              296.23,       // M0 (deg)
              6.2,          // radius (km)
              0.003,        // surface gravity (m/s^2)
              0.00556,      // escape velocity (km/s)
              false, 0.0);

    // ==============================
    // Jupiter's Galilean moons
    // ==============================
    make_moon("Io", "Jupiter",
              8.9319e22, MU_IO, MU_JUPITER,
              421700.0,     // a (km)
              0.0041,       // e
              0.036,        // i (deg) — to Jupiter's equator
              43.98,        // Omega (deg)
              84.13,        // omega (deg)
              342.02,       // M0 (deg)
              1821.6,       // radius (km)
              1.796,        // surface gravity (m/s^2)
              2.558,        // escape velocity (km/s)
              true, 1e-9);  // tenuous SO2 atmosphere

    make_moon("Europa", "Jupiter",
              4.7998e22, MU_EUROPA, MU_JUPITER,
              671100.0,     // a (km)
              0.009,        // e
              0.466,        // i (deg)
              219.11,       // Omega (deg)
              88.97,        // omega (deg)
              171.02,       // M0 (deg)
              1560.8,       // radius (km)
              1.314,        // surface gravity (m/s^2)
              2.025,        // escape velocity (km/s)
              true, 1e-12); // tenuous O2 atmosphere

    make_moon("Ganymede", "Jupiter",
              1.4819e23, MU_GANYMEDE, MU_JUPITER,
              1070400.0,    // a (km)
              0.0013,       // e
              0.177,        // i (deg)
              63.55,        // Omega (deg)
              192.42,       // omega (deg)
              317.54,       // M0 (deg)
              2634.1,       // radius (km)
              1.428,        // surface gravity (m/s^2)
              2.741,        // escape velocity (km/s)
              true, 1e-12); // tenuous O2 atmosphere

    make_moon("Callisto", "Jupiter",
              1.0759e23, MU_CALLISTO, MU_JUPITER,
              1882700.0,    // a (km)
              0.0074,       // e
              0.192,        // i (deg)
              298.85,       // Omega (deg)
              52.64,        // omega (deg)
              181.41,       // M0 (deg)
              2410.3,       // radius (km)
              1.235,        // surface gravity (m/s^2)
              2.440,        // escape velocity (km/s)
              true, 1e-12); // tenuous CO2 atmosphere

    // ==============================
    // Saturn's Titan
    // ==============================
    make_moon("Titan", "Saturn",
              1.3452e23, MU_TITAN, MU_SATURN,
              1221870.0,    // a (km)
              0.0288,       // e
              0.312,        // i (deg) — to Saturn's equator
              169.15,       // Omega (deg)
              185.67,       // omega (deg)
              120.56,       // M0 (deg)
              2574.7,       // radius (km)
              1.352,        // surface gravity (m/s^2)
              2.639,        // escape velocity (km/s)
              true, 1.45);  // thick N2 atmosphere, 1.45 atm
}

} // namespace solar
