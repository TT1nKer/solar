#pragma once

#include <string>
#include <cmath>

namespace solar {

// Propulsion subsystem with uncertainty
struct Propulsion {
    double fuel_mass_kg = 0;
    double thrust_N = 0;
    double isp_s = 0;
    double thrust_uncertainty = 0;  // fractional (0.05 = ±5%)
    double isp_uncertainty = 0;
};

// Vehicle: spacecraft with subsystems
struct Vehicle {
    std::string name;
    double dry_mass_kg = 0;
    Propulsion propulsion;

    // Maximum delta-v from Tsiolkovsky at nominal parameters
    double max_delta_v() const {
        if (propulsion.isp_s <= 0 || dry_mass_kg <= 0) return 0;
        double ve = propulsion.isp_s * 9.80665 / 1000.0; // km/s
        double m0 = dry_mass_kg + propulsion.fuel_mass_kg;
        return ve * std::log(m0 / dry_mass_kg);
    }

    // Max delta-v with perturbed Isp
    double max_delta_v(double isp_factor) const {
        double isp = propulsion.isp_s * isp_factor;
        if (isp <= 0 || dry_mass_kg <= 0) return 0;
        double ve = isp * 9.80665 / 1000.0;
        double m0 = dry_mass_kg + propulsion.fuel_mass_kg;
        return ve * std::log(m0 / dry_mass_kg);
    }
};

} // namespace solar
