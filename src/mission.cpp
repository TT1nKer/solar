#include "solar/mission.h"
#include "solar/ephemeris.h"
#include "solar/kepler.h"
#include "solar/transfer.h"
#include "solar/trajectory.h"
#include "solar/spacecraft.h"
#include "solar/constants.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace solar {

// Compute departure delta-v from parking orbit using vis-viva
static double departure_dv(double v_inf, double mu_body, double r_park) {
    double v_park = std::sqrt(mu_body / r_park);
    double v_hyp = std::sqrt(v_inf * v_inf + 2.0 * mu_body / r_park);
    return v_hyp - v_park;
}

// Compute capture delta-v into circular orbit using vis-viva
static double capture_dv(double v_inf, double mu_body, double r_orbit) {
    double v_circ = std::sqrt(mu_body / r_orbit);
    double v_hyp = std::sqrt(v_inf * v_inf + 2.0 * mu_body / r_orbit);
    return v_hyp - v_circ;
}

// Apply impulsive burn: deduct fuel via Tsiolkovsky
// Returns: {actual_dv, full_burn_achieved}
struct BurnResult { double actual_dv; bool complete; };

static BurnResult apply_burn(double dv_kms, double& fuel_kg, double dry_mass_kg, double isp_s) {
    double ve = isp_s * 9.80665 / 1000.0;
    double m0 = dry_mass_kg + fuel_kg;
    double m1 = m0 * std::exp(-dv_kms / ve);
    double fuel_used = m0 - m1;
    if (fuel_used > fuel_kg + 0.01) {
        double actual_dv = ve * std::log(m0 / dry_mass_kg);
        fuel_kg = 0;
        return {actual_dv, false};
    }
    fuel_kg -= fuel_used;
    return {dv_kms, true};
}

MissionReport execute_mission(const MissionPlan& plan) {
    MissionReport report;
    report.mission_name = plan.name;
    report.dry_mass_kg = plan.dry_mass_kg;
    report.initial_fuel_kg = plan.fuel_mass_kg;
    report.total_dv_kms = 0;
    report.success = true;

    // Max delta-v from Tsiolkovsky
    double ve = plan.isp_s * 9.80665 / 1000.0;
    double m0 = plan.dry_mass_kg + plan.fuel_mass_kg;
    report.max_dv_kms = ve * std::log(m0 / plan.dry_mass_kg);

    auto bodies = load_solar_system();
    double fuel = plan.fuel_mass_kg;
    State sc_state = {{0,0,0},{0,0,0}};
    double current_jd = 0;

    for (const auto& event : plan.events) {
        MissionPhaseResult phase;
        phase.type = event.type;
        phase.jd_start = event.jd;
        phase.description = event.description;
        phase.fuel_before_kg = fuel;
        phase.dv_used_kms = 0;

        switch (event.type) {
        case EventType::Launch: {
            const Body* dep = find_body(bodies, event.departure_body);
            const Body* arr = find_body(bodies, event.arrival_body);
            if (!dep || !arr) {
                report.success = false;
                report.failure_reason = "Unknown body: " + event.departure_body;
                return report;
            }

            State dep_state = compute_position(*dep, event.jd);
            State arr_state = compute_position(*arr, event.arrival_jd);

            double tof = (event.arrival_jd - event.jd) * constants::DAY;
            auto lam = solve_lambert(dep_state.pos, arr_state.pos, tof,
                                     constants::MU_SUN, true);

            Vec3 v_inf_vec = lam.v1 - dep_state.vel;
            double v_inf = v_inf_vec.norm();
            double r_park = dep->properties.radius_km + event.parking_orbit_alt_km;
            double dv = departure_dv(v_inf, dep->mu, r_park);

            auto burn_r = apply_burn(dv, fuel, plan.dry_mass_kg, plan.isp_s);
            double actual_dv = burn_r.actual_dv;
            if (!burn_r.complete) { report.success = false; report.failure_reason = "Insufficient fuel at: " + event.description; }
            phase.dv_used_kms = actual_dv;

            // Set spacecraft state to departure trajectory
            sc_state.pos = dep_state.pos;
            sc_state.vel = lam.v1;

            double c3 = v_inf * v_inf;
            phase.notes = "C3=" + std::to_string(c3).substr(0,5) + " km2/s2";
            phase.jd_end = event.jd;
            current_jd = event.jd;
            break;
        }

        case EventType::Coast: {
            double coast_days = event.coast_end_jd - current_jd;
            phase.notes = std::to_string(static_cast<int>(coast_days)) + " days";

            // If the previous event was a Lambert-based launch, update velocity
            // to the Lambert arrival velocity at the coast end time.
            // Find the last launch event to get arrival info
            for (int k = static_cast<int>(report.phases.size()) - 1; k >= 0; --k) {
                if (report.phases[k].type == EventType::Launch) {
                    // Re-solve Lambert for arrival velocity
                    for (const auto& prev_ev : plan.events) {
                        if (prev_ev.type == EventType::Launch) {
                            const Body* dep = find_body(bodies, prev_ev.departure_body);
                            const Body* arr = find_body(bodies, prev_ev.arrival_body);
                            if (dep && arr) {
                                State dep_s = compute_position(*dep, prev_ev.jd);
                                State arr_s = compute_position(*arr, prev_ev.arrival_jd);
                                double tof = (prev_ev.arrival_jd - prev_ev.jd) * constants::DAY;
                                auto lam = solve_lambert(dep_s.pos, arr_s.pos, tof,
                                                         constants::MU_SUN, true);
                                sc_state.vel = lam.v2; // arrival velocity
                                sc_state.pos = arr_s.pos; // near arrival body
                            }
                            break;
                        }
                    }
                    break;
                }
            }

            phase.jd_end = event.coast_end_jd;
            current_jd = event.coast_end_jd;
            break;
        }

        case EventType::Burn: {
            double dv = event.burn_dv_kms;
            auto burn_r = apply_burn(dv, fuel, plan.dry_mass_kg, plan.isp_s);
            double actual_dv = burn_r.actual_dv;
            if (!burn_r.complete) { report.success = false; report.failure_reason = "Insufficient fuel at: " + event.description; }
            phase.dv_used_kms = actual_dv;

            Vec3 dir = event.burn_direction;
            if (dir.norm() < 0.1) dir = sc_state.vel.normalized();
            sc_state.vel += dir * actual_dv;

            phase.jd_end = event.jd;
            current_jd = event.jd;
            break;
        }

        case EventType::Flyby: {
            const Body* fb = find_body(bodies, event.flyby_body);
            if (!fb) {
                report.success = false;
                report.failure_reason = "Unknown flyby body: " + event.flyby_body;
                return report;
            }

            State fb_state = compute_position(*fb, event.jd);

            FlybyParams fp;
            fp.mu_body = fb->mu;
            fp.radius_body = fb->properties.radius_km;
            fp.v_body = fb_state.vel;
            fp.v_in_helio = sc_state.vel;
            fp.rp = event.flyby_rp_km;

            auto result = compute_flyby(fp);
            if (!result.valid) {
                report.success = false;
                report.failure_reason = "Flyby failed at " + event.flyby_body;
                return report;
            }

            sc_state.vel = result.v_out_helio;
            sc_state.pos = fb_state.pos; // snap to flyby body

            phase.notes = "turn=" + std::to_string(static_cast<int>(
                result.turn_angle * constants::RAD2DEG)) + "deg dv_free=" +
                std::to_string(result.dv_equiv).substr(0,5) + " km/s";
            phase.jd_end = event.jd;
            current_jd = event.jd;
            break;
        }

        case EventType::OrbitInsertion: {
            const Body* tgt = find_body(bodies, event.target_body);
            if (!tgt) {
                report.success = false;
                report.failure_reason = "Unknown target: " + event.target_body;
                return report;
            }

            State tgt_state = compute_position(*tgt, event.jd);
            Vec3 v_inf_vec = sc_state.vel - tgt_state.vel;
            double v_inf = v_inf_vec.norm();
            double r_orbit = tgt->properties.radius_km + event.target_orbit_alt_km;
            double dv = capture_dv(v_inf, tgt->mu, r_orbit);

            auto burn_r = apply_burn(dv, fuel, plan.dry_mass_kg, plan.isp_s);
            double actual_dv = burn_r.actual_dv;
            if (!burn_r.complete) { report.success = false; report.failure_reason = "Insufficient fuel at: " + event.description; }
            phase.dv_used_kms = actual_dv;

            double v_circ = std::sqrt(tgt->mu / r_orbit);
            phase.notes = std::to_string(static_cast<int>(event.target_orbit_alt_km)) +
                          " km orbit, v_circ=" + std::to_string(v_circ).substr(0,4) + " km/s";
            phase.jd_end = event.jd;
            current_jd = event.jd;
            break;
        }
        }

        // No additional check needed — failures detected in event handlers

        phase.fuel_after_kg = fuel;
        phase.state_after = sc_state;
        report.total_dv_kms += phase.dv_used_kms;
        report.phases.push_back(phase);
    }

    report.fuel_remaining_kg = fuel;
    if (!plan.events.empty()) {
        report.total_duration_days = plan.events.back().jd - plan.events.front().jd;
        if (report.total_duration_days < 0) report.total_duration_days = 0;
    }

    return report;
}

void print_mission_report(const MissionReport& report) {
    std::cout << std::string(60, '=') << "\n";
    std::cout << "  MISSION REPORT: " << report.mission_name << "\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << std::fixed;
    std::cout << "  Dry mass: " << std::setprecision(0) << report.dry_mass_kg
              << " kg | Initial fuel: " << report.initial_fuel_kg << " kg\n";
    std::cout << "  Max delta-v: " << std::setprecision(2) << report.max_dv_kms << " km/s\n";
    std::cout << std::string(60, '-') << "\n\n";

    std::cout << std::left
              << std::setw(4) << "#"
              << std::setw(13) << "Date"
              << std::setw(22) << "Event"
              << std::right
              << std::setw(9) << "dv(km/s)"
              << std::setw(10) << "Fuel(kg)"
              << "  Notes\n";
    std::cout << std::string(70, '-') << "\n";

    int phase_num = 1;
    for (const auto& p : report.phases) {
        int y, m, d; double h;
        jd_to_calendar(p.jd_start, y, m, d, h);

        const char* type_str = "";
        switch (p.type) {
            case EventType::Launch: type_str = "Launch"; break;
            case EventType::Coast: type_str = "Coast"; break;
            case EventType::Burn: type_str = "Burn"; break;
            case EventType::Flyby: type_str = "Flyby"; break;
            case EventType::OrbitInsertion: type_str = "Orbit Insertion"; break;
        }

        std::string desc = p.description.empty() ? type_str : p.description;

        std::cout << std::left << std::setw(4) << phase_num
                  << std::right << y << "-" << std::setfill('0') << std::setw(2) << m
                  << "-" << std::setw(2) << d << std::setfill(' ') << "  "
                  << std::left << std::setw(22) << desc
                  << std::right << std::fixed
                  << std::setw(9) << std::setprecision(3) << p.dv_used_kms
                  << std::setw(10) << std::setprecision(0) << p.fuel_after_kg
                  << "  " << p.notes << "\n";
        phase_num++;
    }

    std::cout << std::string(70, '-') << "\n";
    std::cout << std::fixed;
    std::cout << "  Total delta-v:    " << std::setprecision(3) << report.total_dv_kms << " km/s\n";
    std::cout << "  Fuel remaining:   " << std::setprecision(0) << report.fuel_remaining_kg
              << " kg (of " << report.initial_fuel_kg << ")\n";
    std::cout << "  Mission duration: " << std::setprecision(0) << report.total_duration_days << " days\n";
    std::cout << "  Status:           " << (report.success ? "SUCCESS" : "FAILED") << "\n";
    if (!report.success) std::cout << "  Reason:           " << report.failure_reason << "\n";
    std::cout << std::string(60, '=') << "\n";
}

// ============================================================
// Mission Templates
// ============================================================

MissionPlan mars_direct_template() {
    MissionPlan p;
    p.name = "Mars Direct 2026";
    p.dry_mass_kg = 1200; p.fuel_mass_kg = 6200; // tight margin for realistic MC
    p.thrust_N = 44000; p.isp_s = 320;
    p.spacecraft_area_m2 = 12; p.Cd = 2.2; p.Cr = 1.5;

    double jd_dep = calendar_to_jd(2026, 10, 26); // optimal from porkchop
    double jd_arr = calendar_to_jd(2027, 9, 3);  // optimal arrival

    MissionEvent launch;
    launch.type = EventType::Launch;
    launch.jd = jd_dep;
    launch.description = "Earth departure";
    launch.departure_body = "Earth";
    launch.arrival_body = "Mars";
    launch.arrival_jd = jd_arr;
    launch.parking_orbit_alt_km = 300;
    p.events.push_back(launch);

    MissionEvent coast;
    coast.type = EventType::Coast;
    coast.jd = jd_dep;
    coast.description = "Interplanetary cruise";
    coast.coast_end_jd = jd_arr;
    p.events.push_back(coast);

    MissionEvent moi;
    moi.type = EventType::OrbitInsertion;
    moi.jd = jd_arr;
    moi.description = "Mars orbit insertion";
    moi.target_body = "Mars";
    moi.target_orbit_alt_km = 400;
    p.events.push_back(moi);

    return p;
}

// LIMITATION: This template uses HARDCODED Δv values (3.13 km/s TLI, 0.82 km/s LOI).
// It does NOT compute the actual Earth-Moon trajectory because:
//   1. Earth-Moon is NOT a heliocentric Lambert problem (Moon is in Earth's SOI)
//   2. Proper computation requires patched-conic Earth-Moon transfer
//      (Earth-centered hyperbolic departure -> Moon SOI -> Moon-centered capture)
//   3. Or three-body dynamics in the Earth-Moon CR3BP
// The Δv values are correct for typical low lunar orbit insertion from LEO,
// taken from textbook references. This is a Δv accounting demonstration only.
MissionPlan lunar_gateway_template() {
    MissionPlan p;
    p.name = "Lunar Gateway (Δv-accounting demo)";
    p.dry_mass_kg = 2000; p.fuel_mass_kg = 6000;
    p.thrust_N = 25000; p.isp_s = 320;
    p.spacecraft_area_m2 = 15; p.Cd = 2.2; p.Cr = 1.5;

    double jd_dep = calendar_to_jd(2026, 6, 15);
    double jd_arr = jd_dep + 4.5; // ~4.5 days transit (standard lunar transfer)

    // TLI: hardcoded standard value from textbook (NOT computed from ephemeris)
    MissionEvent tli;
    tli.type = EventType::Burn;
    tli.jd = jd_dep;
    tli.description = "Trans-Lunar Injection";
    tli.burn_dv_kms = 3.13; // textbook TLI from 200 km LEO
    tli.burn_direction = {0, 0, 0};
    p.events.push_back(tli);

    MissionEvent coast;
    coast.type = EventType::Coast;
    coast.jd = jd_dep;
    coast.description = "Lunar transit (no propagation)";
    coast.coast_end_jd = jd_arr;
    p.events.push_back(coast);

    // LOI: hardcoded standard value (NOT computed from ephemeris)
    MissionEvent loi;
    loi.type = EventType::Burn;
    loi.jd = jd_arr;
    loi.description = "Lunar orbit insertion";
    loi.burn_dv_kms = 0.82; // textbook LOI value
    loi.burn_direction = {0, 0, 0};
    p.events.push_back(loi);

    return p;
}

MissionPlan grand_tour_template() {
    MissionPlan p;
    p.name = "Grand Tour (Voyager-like)";
    p.dry_mass_kg = 475; p.fuel_mass_kg = 425; // includes Centaur upper stage contribution
    p.thrust_N = 440; p.isp_s = 230;
    p.spacecraft_area_m2 = 10; p.Cd = 2.2; p.Cr = 1.5;

    double jd_dep = calendar_to_jd(1977, 9, 5);
    double jd_jup = calendar_to_jd(1979, 7, 9);   // Jupiter flyby
    double jd_sat = calendar_to_jd(1981, 8, 25);   // Saturn flyby

    MissionEvent launch;
    launch.type = EventType::Launch;
    launch.jd = jd_dep;
    launch.description = "Earth departure";
    launch.departure_body = "Earth";
    launch.arrival_body = "Jupiter";
    launch.arrival_jd = jd_jup;
    launch.parking_orbit_alt_km = 300;
    p.events.push_back(launch);

    MissionEvent coast1;
    coast1.type = EventType::Coast;
    coast1.jd = jd_dep;
    coast1.description = "Earth-Jupiter cruise";
    coast1.coast_end_jd = jd_jup;
    p.events.push_back(coast1);

    MissionEvent jup_flyby;
    jup_flyby.type = EventType::Flyby;
    jup_flyby.jd = jd_jup;
    jup_flyby.description = "Jupiter flyby";
    jup_flyby.flyby_body = "Jupiter";
    jup_flyby.flyby_rp_km = 350000; // ~5 Jupiter radii
    p.events.push_back(jup_flyby);

    MissionEvent coast2;
    coast2.type = EventType::Coast;
    coast2.jd = jd_jup;
    coast2.description = "Jupiter-Saturn cruise";
    coast2.coast_end_jd = jd_sat;
    p.events.push_back(coast2);

    MissionEvent sat_flyby;
    sat_flyby.type = EventType::Flyby;
    sat_flyby.jd = jd_sat;
    sat_flyby.description = "Saturn flyby";
    sat_flyby.flyby_body = "Saturn";
    sat_flyby.flyby_rp_km = 180000; // ~3 Saturn radii
    p.events.push_back(sat_flyby);

    return p;
}

} // namespace solar
