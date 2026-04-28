#pragma once

#include "solar/body.h"

namespace solar {

struct TransferResult {
    double dv1;           // departure delta-v (km/s)
    double dv2;           // arrival delta-v (km/s)
    double dv_total;      // total delta-v (km/s)
    double tof;           // time of flight (seconds)
    double a_transfer;    // semi-major axis of transfer orbit (km)
    double e_transfer;    // eccentricity of transfer orbit
};

// Hohmann transfer between two circular coplanar orbits
TransferResult hohmann_transfer(double r1, double r2, double mu);

// Bi-elliptic transfer (can be cheaper than Hohmann when r2/r1 > 11.94)
TransferResult bielliptic_transfer(double r1, double r2, double r_intermediate, double mu);

// Lambert problem result
struct LambertResult {
    Vec3 v1;          // departure velocity vector
    Vec3 v2;          // arrival velocity vector
    int iterations;   // solver iterations used
    bool converged;
};

// Solve Lambert's problem: find orbit connecting r1 to r2 in time tof
// Uses universal variable formulation
LambertResult solve_lambert(const Vec3& r1, const Vec3& r2, double tof, double mu,
                            bool prograde = true);

// Plan a transfer between two bodies at given dates
// Returns delta-v requirements using Lambert solution
TransferResult plan_transfer(const Body& departure, const Body& arrival,
                             double jd_depart, double jd_arrive, double mu_central);

} // namespace solar
