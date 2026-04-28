#pragma once

#include <string>
#include <cmath>

namespace solar {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    double norm() const { return std::sqrt(x * x + y * y + z * z); }
    double norm_sq() const { return x * x + y * y + z * z; }
    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    Vec3 normalized() const { double n = norm(); return {x / n, y / n, z / n}; }
};

inline Vec3 operator*(double s, const Vec3& v) { return {s * v.x, s * v.y, s * v.z}; }

// Cartesian state vector
struct State {
    Vec3 pos;   // km
    Vec3 vel;   // km/s
};

// Classical Keplerian orbital elements
struct OrbitalElements {
    double a;       // semi-major axis (km)
    double e;       // eccentricity
    double i;       // inclination (rad)
    double omega;   // argument of periapsis (rad)
    double Omega;   // longitude of ascending node (rad)
    double M0;      // mean anomaly at epoch (rad)
    double epoch;   // reference epoch (Julian date)
};

enum class BodyType { Star, Planet, Moon, Spacecraft };

struct PhysicalProperties {
    double radius_km = 0.0;          // mean volumetric radius
    double surface_gravity = 0.0;    // m/s^2
    double escape_velocity = 0.0;    // km/s
    bool has_atmosphere = false;
    double atm_surface_pressure = 0.0; // atm (0 = none)
};

// A celestial body
struct Body {
    std::string name;
    double mass;                // kg
    double mu;                  // gravitational parameter GM (km^3/s^2)
    OrbitalElements elements;   // osculating elements at epoch
    State state;                // current cartesian state

    // Phase 2: hierarchy and physical properties
    BodyType type = BodyType::Planet;
    int parent_index = -1;      // index into bodies vector; -1 = root (Sun)
    double mu_parent = 0.0;     // GM of parent body (for Keplerian propagation)
    PhysicalProperties properties;
};

} // namespace solar
