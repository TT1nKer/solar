#pragma once

#include "solar/body.h"
#include <string>

namespace solar {

// 3x3 rotation matrix (row-major)
struct Mat3 {
    double m[3][3] = {};

    Vec3 apply(const Vec3& v) const {
        return {
            m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
            m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
            m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
        };
    }

    Mat3 operator*(const Mat3& B) const {
        Mat3 R = {};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    R.m[i][j] += m[i][k] * B.m[k][j];
        return R;
    }

    Mat3 transpose() const {
        Mat3 R;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                R.m[i][j] = m[j][i];
        return R;
    }

    static Mat3 identity() {
        Mat3 R = {};
        R.m[0][0] = R.m[1][1] = R.m[2][2] = 1.0;
        return R;
    }

    static Mat3 rotation_x(double a);
    static Mat3 rotation_y(double a);
    static Mat3 rotation_z(double a);
};

// Reference frames
enum class Frame {
    EclipticJ2000,   // Current default: heliocentric ecliptic J2000
    ICRF,            // Equatorial J2000 (= ICRF). DE440 native frame.
    BodyFixed,       // IAU rotating body-fixed frame (requires body_id + time)
};

// Parse frame name string (case-insensitive)
Frame parse_frame(const std::string& name);
const char* frame_name(Frame f);

// Constant rotation matrices
const Mat3& ecliptic_to_icrf();
const Mat3& icrf_to_ecliptic();

// Transform a state between frames
// For BodyFixed, body_id is required (DE convention: 3=Earth, 4=Mars, 10=Moon)
State transform_state(const State& s, Frame from, Frame to,
                      double jd_tdb, int body_id = -1);

} // namespace solar
