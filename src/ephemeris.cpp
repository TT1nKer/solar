#include "solar/ephemeris.h"
#include "solar/moons.h"
#include "solar/kepler.h"
#include "solar/jpl_de.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace solar {

// JPL mean orbital elements at J2000.0 (ecliptic and mean equinox of J2000)
// Source: Standish (1992), JPL Solar System Dynamics
// Units: a in AU (converted to km), angles in degrees (converted to radians)
std::vector<Body> load_solar_system() {
    using namespace constants;

    std::vector<Body> bodies;

    // Sun (at origin in heliocentric frame)
    {
        Body sun;
        sun.name = "Sun";
        sun.mass = 1.989e30;
        sun.mu = MU_SUN;
        sun.elements = {0, 0, 0, 0, 0, 0, J2000};
        sun.state = {{0, 0, 0}, {0, 0, 0}};
        sun.type = BodyType::Star;
        sun.parent_index = -1;
        sun.properties = {695700.0, 274.0, 617.7, false, 0.0};
        bodies.push_back(sun);
    }

    // Helper: create planet from JPL elements
    // a_au: semi-major axis (AU)
    // e: eccentricity
    // i_deg: inclination (deg)
    // Omega_deg: lon ascending node (deg)
    // omega_bar_deg: longitude of perihelion (deg)
    // L_deg: mean longitude (deg)
    auto make_planet = [&](const std::string& name, double mass_kg, double mu_val,
                           double a_au, double e, double i_deg,
                           double Omega_deg, double omega_bar_deg, double L_deg,
                           PhysicalProperties props = {}) {
        Body b;
        b.name = name;
        b.mass = mass_kg;
        b.mu = mu_val;
        b.type = BodyType::Planet;
        b.parent_index = 0; // Sun is always index 0
        b.mu_parent = MU_SUN;
        b.properties = props;

        double omega_bar = omega_bar_deg * DEG2RAD;
        double Omega = Omega_deg * DEG2RAD;
        double omega = omega_bar - Omega; // argument of periapsis
        double L = L_deg * DEG2RAD;
        double M0 = L - omega_bar; // mean anomaly

        b.elements.a = a_au * AU;
        b.elements.e = e;
        b.elements.i = i_deg * DEG2RAD;
        b.elements.omega = omega;
        b.elements.Omega = Omega;
        b.elements.M0 = M0;
        b.elements.epoch = J2000;

        // Compute initial state at J2000
        b.state = elements_to_state_at(b.elements, MU_SUN, J2000);

        bodies.push_back(b);
    };

    //                                                                          radius   g      v_esc  atmo   P_atm
    // Mercury
    make_planet("Mercury", 3.301e23, MU_MERCURY,
                0.38709927, 0.20563593, 7.00497902,
                48.33076593, 77.45779628, 252.25032350,
                {2439.7, 3.70, 4.25, false, 0.0});

    // Venus
    make_planet("Venus", 4.867e24, MU_VENUS,
                0.72333566, 0.00677672, 3.39467605,
                76.67984255, 131.60246718, 181.97909950,
                {6051.8, 8.87, 10.36, true, 92.0});

    // Earth
    make_planet("Earth", 5.972e24, MU_EARTH,
                1.00000261, 0.01671123, 0.00001531,
                -11.26064, 102.93768193, 100.46457166,
                {6371.0, 9.81, 11.186, true, 1.0});

    // Mars
    make_planet("Mars", 6.417e23, MU_MARS,
                1.52371034, 0.09339410, 1.84969142,
                49.55953891, 336.05637041, 355.44656895,
                {3389.5, 3.72, 5.03, true, 0.006});

    // Jupiter
    make_planet("Jupiter", 1.898e27, MU_JUPITER,
                5.20288700, 0.04838624, 1.30439695,
                100.47390909, 14.72847983, 34.39644051,
                {69911.0, 24.79, 59.5, true, 0.0});  // no solid surface

    // Saturn
    make_planet("Saturn", 5.683e26, MU_SATURN,
                9.53667594, 0.05386179, 2.48599187,
                113.66242448, 92.59887831, 49.95424423,
                {58232.0, 10.44, 35.5, true, 0.0});

    // Uranus
    make_planet("Uranus", 8.681e25, MU_URANUS,
                19.18916464, 0.04725744, 0.77263783,
                74.01692503, 170.95427630, 313.23810451,
                {25362.0, 8.87, 21.3, true, 0.0});

    // Neptune
    make_planet("Neptune", 1.024e26, MU_NEPTUNE,
                30.06992276, 0.00859048, 1.77004347,
                131.78422574, 44.96476227, 304.88003433,
                {24622.0, 11.15, 23.5, true, 0.0});

    // Load moons (must happen after planets are in the vector)
    load_moons(bodies);

    return bodies;
}

const Body* find_body(const std::vector<Body>& bodies, const std::string& name) {
    std::string lower_name = name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    for (const auto& b : bodies) {
        std::string lower_b = b.name;
        std::transform(lower_b.begin(), lower_b.end(), lower_b.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (lower_b == lower_name) return &b;
    }
    return nullptr;
}

State compute_position(const Body& body, double jd) {
    if (body.type == BodyType::Star || body.parent_index < 0) {
        return {{0, 0, 0}, {0, 0, 0}};
    }

    // Try JPL DE ephemeris for planets and Moon
    auto* de = global_de_ephemeris();
    if (de && de->covers(jd)) {
        int de_id = body_name_to_de_id(body.name);
        if (de_id > 0) {
            State s;
            if (de->get_state(de_id, jd, s)) {
                return s;
            }
        }
    }

    // Fallback: Keplerian propagation
    double mu = (body.mu_parent > 0.0) ? body.mu_parent : constants::MU_SUN;
    return elements_to_state_at(body.elements, mu, jd);
}

State compute_heliocentric(const Body& body, const std::vector<Body>& bodies, double jd) {
    State s = compute_position(body, jd);
    int pi = body.parent_index;
    while (pi > 0) { // walk up until Sun (index 0)
        State parent_s = compute_position(bodies[pi], jd);
        s.pos = s.pos + parent_s.pos;
        s.vel = s.vel + parent_s.vel;
        pi = bodies[pi].parent_index;
    }
    return s;
}

double calendar_to_jd(int year, int month, int day, double hour) {
    // Algorithm from Meeus, "Astronomical Algorithms"
    if (month <= 2) {
        year -= 1;
        month += 12;
    }
    int A = year / 100;
    int B = 2 - A + A / 4;

    double JD = std::floor(365.25 * (year + 4716))
              + std::floor(30.6001 * (month + 1))
              + day + hour / 24.0 + B - 1524.5;
    return JD;
}

void jd_to_calendar(double jd, int& year, int& month, int& day, double& hour) {
    double z_d = std::floor(jd + 0.5);
    double f = (jd + 0.5) - z_d;
    int Z = static_cast<int>(z_d);

    int A;
    if (Z < 2299161) {
        A = Z;
    } else {
        int alpha = static_cast<int>((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - alpha / 4;
    }

    int B = A + 1524;
    int C = static_cast<int>((B - 122.1) / 365.25);
    int D = static_cast<int>(365.25 * C);
    int E = static_cast<int>((B - D) / 30.6001);

    day = B - D - static_cast<int>(30.6001 * E);
    month = (E < 14) ? E - 1 : E - 13;
    year = (month > 2) ? C - 4716 : C - 4715;
    hour = f * 24.0;
}

double parse_date(const std::string& date_str) {
    int year, month, day;
    char sep1, sep2;
    std::istringstream iss(date_str);
    iss >> year >> sep1 >> month >> sep2 >> day;
    if (iss.fail()) {
        throw std::runtime_error("Invalid date format. Expected YYYY-MM-DD");
    }
    return calendar_to_jd(year, month, day);
}

} // namespace solar
