#pragma once

#include "solar/body.h"
#include <vector>
#include <string>

namespace solar {

// Load the solar system: Sun + 8 planets + moons with J2000 orbital elements
std::vector<Body> load_solar_system();

// Find a body by name (case-insensitive)
const Body* find_body(const std::vector<Body>& bodies, const std::string& name);

// Compute state of a body relative to its parent at a given Julian date
// For planets: heliocentric. For moons: relative to parent planet.
State compute_position(const Body& body, double jd);

// Compute heliocentric state by walking up the parent chain
// Works for any body (planets, moons, spacecraft)
State compute_heliocentric(const Body& body, const std::vector<Body>& bodies, double jd);

// Convert calendar date to Julian date
double calendar_to_jd(int year, int month, int day, double hour = 0.0);

// Convert Julian date to calendar date
void jd_to_calendar(double jd, int& year, int& month, int& day, double& hour);

// Parse date string "YYYY-MM-DD" to Julian date
double parse_date(const std::string& date_str);

} // namespace solar
