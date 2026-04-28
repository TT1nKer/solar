#include "solar/iau_rotation.h"
#include "solar/constants.h"
#include <cmath>
#include <stdexcept>

namespace solar {

// IAU 2015 rotation models
// Reference: Archinal et al. (2018), "Report of the IAU Working Group on
// Cartographic Coordinates and Rotational Elements: 2015"

static const IAURotation IAU_EARTH = {
    0.0, -0.641,           // RA of pole (degrees + degrees/century)
    90.0, -0.557,          // Dec of pole
    190.147, 360.9856235   // Prime meridian (degrees + degrees/day)
};

static const IAURotation IAU_MARS = {
    317.68143, -0.1061,
    52.88650, -0.0609,
    176.630, 350.89198226
};

static const IAURotation IAU_MOON = {
    // Simplified (without periodic terms E1-E13)
    269.9949, 0.0031,
    66.5392, 0.0130,
    38.3213, 13.17635815
};

// Body ID to IAU model
static const IAURotation IAU_DEFAULT = {0, 0, 90, 0, 0, 0};

const IAURotation& iau_rotation_params(int body_id) {
    switch (body_id) {
        case 3:  return IAU_EARTH;
        case 4:  return IAU_MARS;
        case 10: return IAU_MOON;
    }
    return IAU_DEFAULT;
}

// Compute ICRF -> body-fixed rotation matrix
// R = Rz(-W) * Rx(-(90-delta)) * Rz(-(90+alpha))
// omega is the angular velocity of the body-fixed frame in ICRF (rad/s)
Mat3 iau_body_fixed_matrix(int body_id, double jd_tdb, Vec3& omega_icrf) {
    using namespace constants;

    const IAURotation& m = iau_rotation_params(body_id);

    double d = jd_tdb - J2000;         // days since J2000
    double T = d / 36525.0;            // Julian centuries since J2000

    // Pole direction (degrees)
    double alpha = m.alpha0 + m.alpha1 * T;
    double delta = m.delta0 + m.delta1 * T;

    // Prime meridian angle (degrees)
    double W = m.W0 + m.W_dot * d;

    // Convert to radians
    double alpha_rad = alpha * DEG2RAD;
    double delta_rad = delta * DEG2RAD;
    double W_rad = W * DEG2RAD;

    // Rotation matrix: ICRF -> body-fixed
    // R = Rz(-W) * Rx(delta - 90°) * Rz(alpha + 90°)
    // This is equivalent to: first align pole, then rotate to prime meridian
    Mat3 R = Mat3::rotation_z(-W_rad)
           * Mat3::rotation_x(delta_rad - PI / 2.0)
           * Mat3::rotation_z(alpha_rad + PI / 2.0);

    // Angular velocity of the body (in ICRF frame)
    // The body's pole direction in ICRF:
    double cos_a = std::cos(alpha_rad), sin_a = std::sin(alpha_rad);
    double cos_d = std::cos(delta_rad), sin_d = std::sin(delta_rad);

    // Pole unit vector in ICRF
    Vec3 pole = {cos_a * cos_d, sin_a * cos_d, sin_d};

    // W_dot in rad/s
    double W_dot_rad_per_sec = m.W_dot * DEG2RAD / DAY;

    // omega = W_dot * pole_direction
    omega_icrf = pole * W_dot_rad_per_sec;

    return R;
}

Mat3 iau_body_fixed_matrix(int body_id, double jd_tdb) {
    Vec3 omega;
    return iau_body_fixed_matrix(body_id, jd_tdb, omega);
}

LatLonAlt cartesian_to_geodetic(const Vec3& pos, double R_eq) {
    double r = pos.norm();
    double lat = std::asin(pos.z / r) * constants::RAD2DEG;
    double lon = std::atan2(pos.y, pos.x) * constants::RAD2DEG;
    double alt = r - R_eq;
    return {lat, lon, alt};
}

} // namespace solar
