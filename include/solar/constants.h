#pragma once

namespace solar {
namespace constants {

constexpr double G       = 6.67430e-20;        // km^3 / (kg * s^2)
constexpr double AU      = 1.495978707e8;       // km
constexpr double PI      = 3.14159265358979323846;
constexpr double TWO_PI  = 2.0 * PI;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double DAY     = 86400.0;             // seconds
constexpr double YEAR    = 365.25 * DAY;        // Julian year in seconds

// Gravitational parameters (km^3/s^2) — known more precisely than G*M
constexpr double MU_SUN     = 1.32712440018e11;
constexpr double MU_MERCURY = 2.2032e4;
constexpr double MU_VENUS   = 3.24859e5;
constexpr double MU_EARTH   = 3.986004418e5;
constexpr double MU_MARS    = 4.282837e4;
constexpr double MU_JUPITER = 1.26686534e8;
constexpr double MU_SATURN  = 3.7931187e7;
constexpr double MU_URANUS  = 5.793939e6;
constexpr double MU_NEPTUNE = 6.836529e6;

// Moon gravitational parameters (km^3/s^2)
constexpr double MU_MOON     = 4.9028695e3;   // Luna
constexpr double MU_IO       = 5.959916e3;
constexpr double MU_EUROPA   = 3.202739e3;
constexpr double MU_GANYMEDE = 9.887834e3;
constexpr double MU_CALLISTO = 7.179289e3;
constexpr double MU_TITAN    = 8.97814e3;
constexpr double MU_PHOBOS   = 7.087546e-4;
constexpr double MU_DEIMOS   = 9.8e-5;

// Speed of light (km/s)
constexpr double C_LIGHT = 299792.458;

// Mean obliquity of ecliptic at J2000 (IAU 2000)
constexpr double OBLIQUITY_J2000 = 23.4392911111 * DEG2RAD;

// Julian date of J2000.0 epoch (2000-01-01 12:00 TT)
constexpr double J2000 = 2451545.0;

} // namespace constants
} // namespace solar
