#pragma once

#include "solar/body.h"
#include <vector>
#include <string>

namespace solar {

enum class EventType {
    Launch,
    Coast,
    Burn,
    Flyby,
    OrbitInsertion,
};

struct MissionEvent {
    EventType type;
    double jd;
    std::string description;

    // Launch
    std::string departure_body;
    std::string arrival_body;     // for Lambert computation
    double arrival_jd;            // arrival date for Lambert
    double parking_orbit_alt_km;

    // Coast
    double coast_end_jd;

    // Burn (impulsive)
    double burn_dv_kms;
    Vec3 burn_direction;          // {0,0,0} = prograde

    // Flyby
    std::string flyby_body;
    double flyby_rp_km;

    // Orbit Insertion
    std::string target_body;
    double target_orbit_alt_km;
};

struct MissionPhaseResult {
    EventType type;
    double jd_start, jd_end;
    std::string description;
    double dv_used_kms;
    double fuel_before_kg;
    double fuel_after_kg;
    State state_after;
    std::string notes;
};

struct MissionReport {
    std::string mission_name;
    double dry_mass_kg;
    double initial_fuel_kg;
    double max_dv_kms;
    std::vector<MissionPhaseResult> phases;
    double total_dv_kms;
    double fuel_remaining_kg;
    double total_duration_days;
    bool success;
    std::string failure_reason;
};

struct MissionPlan {
    std::string name;
    double dry_mass_kg;
    double fuel_mass_kg;
    double thrust_N;
    double isp_s;
    double spacecraft_area_m2;
    double Cd, Cr;
    std::vector<MissionEvent> events;
};

// Execute a mission plan and return the report
MissionReport execute_mission(const MissionPlan& plan);

// Print formatted mission report
void print_mission_report(const MissionReport& report);

// Pre-built mission templates
MissionPlan mars_direct_template();
MissionPlan lunar_gateway_template();
MissionPlan grand_tour_template();

} // namespace solar
