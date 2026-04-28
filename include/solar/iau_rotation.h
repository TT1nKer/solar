#pragma once

#include "solar/frame.h"

namespace solar {

// IAU rotation model parameters for a body
struct IAURotation {
    double alpha0, alpha1;   // RA of pole: alpha0 + alpha1*T (degrees)
    double delta0, delta1;   // Dec of pole: delta0 + delta1*T (degrees)
    double W0, W_dot;        // Prime meridian: W0 + W_dot*d (degrees)
};

// Get IAU rotation parameters for a body
// body_id: 3=Earth, 4=Mars, 10=Moon, etc.
const IAURotation& iau_rotation_params(int body_id);

// Compute ICRF-to-body-fixed rotation matrix at given TDB epoch
// Also returns angular velocity vector in ICRF (rad/s) for velocity transforms
Mat3 iau_body_fixed_matrix(int body_id, double jd_tdb, Vec3& omega_icrf);

// Convenience: just the matrix (no omega)
Mat3 iau_body_fixed_matrix(int body_id, double jd_tdb);

// Convert body-fixed Cartesian to latitude/longitude/altitude
// lat: degrees (-90 to 90), lon: degrees (-180 to 180), alt: km above sphere
struct LatLonAlt {
    double lat_deg;
    double lon_deg;
    double alt_km;
};

LatLonAlt cartesian_to_geodetic(const Vec3& body_fixed_pos, double R_eq);

} // namespace solar
