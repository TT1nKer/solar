#include "solar/frame.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>
#include <cctype>
#include <stdexcept>

// Forward declaration — implemented in iau_rotation.cpp
namespace solar {
    Mat3 iau_body_fixed_matrix(int body_id, double jd_tdb, Vec3& omega_icrf);
}

namespace solar {

Mat3 Mat3::rotation_x(double a) {
    double c = std::cos(a), s = std::sin(a);
    Mat3 R = {};
    R.m[0][0] = 1.0;
    R.m[1][1] = c;  R.m[1][2] = -s;
    R.m[2][1] = s;  R.m[2][2] = c;
    return R;
}

Mat3 Mat3::rotation_y(double a) {
    double c = std::cos(a), s = std::sin(a);
    Mat3 R = {};
    R.m[0][0] = c;  R.m[0][2] = s;
    R.m[1][1] = 1.0;
    R.m[2][0] = -s; R.m[2][2] = c;
    return R;
}

Mat3 Mat3::rotation_z(double a) {
    double c = std::cos(a), s = std::sin(a);
    Mat3 R = {};
    R.m[0][0] = c;  R.m[0][1] = -s;
    R.m[1][0] = s;  R.m[1][1] = c;
    R.m[2][2] = 1.0;
    return R;
}

// Constant ecliptic ↔ ICRF rotation (around X by obliquity)
static Mat3 make_ecl_to_icrf() {
    return Mat3::rotation_x(-constants::OBLIQUITY_J2000);
}

static Mat3 make_icrf_to_ecl() {
    return Mat3::rotation_x(constants::OBLIQUITY_J2000);
}

const Mat3& ecliptic_to_icrf() {
    static const Mat3 R = make_ecl_to_icrf();
    return R;
}

const Mat3& icrf_to_ecliptic() {
    static const Mat3 R = make_icrf_to_ecl();
    return R;
}

Frame parse_frame(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (lower == "ecliptic" || lower == "eclipticj2000") return Frame::EclipticJ2000;
    if (lower == "icrf" || lower == "equatorial" || lower == "j2000") return Frame::ICRF;
    if (lower == "bodyfixed" || lower == "body_fixed" || lower == "fixed") return Frame::BodyFixed;
    throw std::runtime_error("Unknown frame: " + name);
}

const char* frame_name(Frame f) {
    switch (f) {
        case Frame::EclipticJ2000: return "EclipticJ2000";
        case Frame::ICRF: return "ICRF";
        case Frame::BodyFixed: return "BodyFixed";
    }
    return "?";
}

// Transform state between frames. All conversions go through ICRF as hub.
State transform_state(const State& s, Frame from, Frame to,
                      double jd_tdb, int body_id) {
    if (from == to) return s;

    // Step 1: Convert to ICRF
    State icrf = s;
    if (from == Frame::EclipticJ2000) {
        const Mat3& R = ecliptic_to_icrf();
        icrf.pos = R.apply(s.pos);
        icrf.vel = R.apply(s.vel);
    } else if (from == Frame::BodyFixed) {
        Vec3 omega;
        Mat3 R = iau_body_fixed_matrix(body_id, jd_tdb, omega);
        Mat3 Rt = R.transpose(); // body-fixed to ICRF
        icrf.pos = Rt.apply(s.pos);
        // Velocity: undo rotating frame correction
        // v_inertial = Rt * (v_body + omega x r_body)
        // But omega is in ICRF, so: v_inertial = Rt * v_body + omega x (Rt * r_body)
        Vec3 r_icrf = icrf.pos;
        icrf.vel = Rt.apply(s.vel) + omega.cross(r_icrf);
    }

    // Step 2: Convert from ICRF to target
    if (to == Frame::ICRF) return icrf;

    State result = icrf;
    if (to == Frame::EclipticJ2000) {
        const Mat3& R = icrf_to_ecliptic();
        result.pos = R.apply(icrf.pos);
        result.vel = R.apply(icrf.vel);
    } else if (to == Frame::BodyFixed) {
        Vec3 omega;
        Mat3 R = iau_body_fixed_matrix(body_id, jd_tdb, omega);
        result.pos = R.apply(icrf.pos);
        // v_body = R * v_inertial - (R * omega) x (R * r) = R * (v_inertial - omega x r)
        Vec3 omega_cross_r = omega.cross(icrf.pos);
        result.vel = R.apply(icrf.vel - omega_cross_r);
    }

    return result;
}

} // namespace solar
