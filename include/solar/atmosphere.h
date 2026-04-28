#pragma once

#include <string>

namespace solar {

// Exponential atmosphere density profile
struct AtmosphereProfile {
    std::string body_name;
    double surface_density;   // kg/m^3 at surface
    double scale_height;      // km
    double surface_radius;    // km (body radius)
    double max_altitude;      // km (above this, density = 0)
};

// Built-in profiles
const AtmosphereProfile& earth_atmosphere();
const AtmosphereProfile& mars_atmosphere();
const AtmosphereProfile& venus_atmosphere();
const AtmosphereProfile& titan_atmosphere();

// Get profile by body name (returns nullptr if no atmosphere model)
const AtmosphereProfile* atmosphere_for_body(const std::string& name);

// Compute density at altitude above surface (km). Returns 0 if out of range.
double atmospheric_density(const AtmosphereProfile& atm, double altitude_km);

} // namespace solar
