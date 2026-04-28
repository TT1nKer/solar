#pragma once

#include "solar/body.h"
#include "solar/transfer.h"
#include <vector>
#include <string>

namespace solar {

// ============================================================
// Porkchop Plot
// ============================================================

struct PorkchopParams {
    std::string departure_body;
    std::string arrival_body;
    double jd_depart_start, jd_depart_end;
    double jd_arrive_start, jd_arrive_end;
    int n_depart = 100;
    int n_arrive = 100;
    double mu_central = 0;  // 0 = MU_SUN
    bool prograde = true;
};

struct PorkchopCell {
    double jd_depart, jd_arrive;
    double dv1, dv2, dv_total;
    double tof_days;
    bool valid;
};

struct PorkchopGrid {
    std::vector<double> depart_dates;
    std::vector<double> arrive_dates;
    std::vector<std::vector<PorkchopCell>> cells; // [depart][arrive]
};

struct LaunchWindow {
    double jd_depart, jd_arrive;
    double dv_total, dv1, dv2;
    double tof_days;
};

PorkchopGrid compute_porkchop(const PorkchopParams& params,
                              const std::vector<Body>& bodies);

LaunchWindow find_optimal_window(const PorkchopGrid& grid);

std::vector<LaunchWindow> find_windows_per_synodic(
    const PorkchopGrid& grid, double synodic_period_days);

// ============================================================
// Gravity Assist
// ============================================================

struct FlybyParams {
    double mu_body;
    double radius_body;   // km
    Vec3 v_body;          // heliocentric velocity of flyby body (km/s)
    Vec3 v_in_helio;      // incoming spacecraft velocity (heliocentric, km/s)
    double rp;            // periapsis distance from body center (km)
};

struct FlybyResult {
    Vec3 v_out_helio;
    Vec3 v_inf_in, v_inf_out;
    double v_inf_mag;
    double turn_angle;    // rad
    double dv_equiv;      // free delta-v from flyby (km/s)
    double e_hyp;         // hyperbolic eccentricity
    double altitude_km;
    bool valid;
};

FlybyResult compute_flyby(const FlybyParams& params);

// ============================================================
// Multi-flyby Trajectory
// ============================================================

struct FlybyLeg {
    std::string body_name;
    double jd;
};

struct MultiFlybyParams {
    std::string departure_body;
    std::string arrival_body;
    double jd_depart, jd_arrive;
    std::vector<FlybyLeg> flybys;
    double mu_central = 0; // 0 = MU_SUN
};

struct MultiFlybyLegResult {
    std::string from_body, to_body;
    double jd_depart, jd_arrive;
    double tof_days;
    double dv;
};

struct MultiFlybyResult {
    std::vector<MultiFlybyLegResult> legs;
    double dv_departure, dv_arrival, dv_total;
    double total_tof_days;
    bool valid;
};

MultiFlybyResult evaluate_multi_flyby(const MultiFlybyParams& params,
                                      const std::vector<Body>& bodies);

// Synodic period between two bodies (days)
double synodic_period(const Body& a, const Body& b, double mu_central);

} // namespace solar
