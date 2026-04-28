#include "solar/atmosphere.h"
#include <cmath>
#include <algorithm>
#include <cctype>

namespace solar {

static const AtmosphereProfile ATM_EARTH = {
    "Earth", 1.225, 8.5, 6371.0, 1000.0
};

static const AtmosphereProfile ATM_MARS = {
    "Mars", 0.020, 11.1, 3389.5, 300.0
};

static const AtmosphereProfile ATM_VENUS = {
    "Venus", 65.0, 15.9, 6051.8, 250.0
};

static const AtmosphereProfile ATM_TITAN = {
    "Titan", 5.3, 21.0, 2574.7, 600.0
};

const AtmosphereProfile& earth_atmosphere() { return ATM_EARTH; }
const AtmosphereProfile& mars_atmosphere()  { return ATM_MARS; }
const AtmosphereProfile& venus_atmosphere() { return ATM_VENUS; }
const AtmosphereProfile& titan_atmosphere() { return ATM_TITAN; }

const AtmosphereProfile* atmosphere_for_body(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (lower == "earth") return &ATM_EARTH;
    if (lower == "mars")  return &ATM_MARS;
    if (lower == "venus") return &ATM_VENUS;
    if (lower == "titan") return &ATM_TITAN;
    return nullptr;
}

double atmospheric_density(const AtmosphereProfile& atm, double altitude_km) {
    if (altitude_km < 0.0 || altitude_km > atm.max_altitude) return 0.0;
    return atm.surface_density * std::exp(-altitude_km / atm.scale_height);
}

} // namespace solar
