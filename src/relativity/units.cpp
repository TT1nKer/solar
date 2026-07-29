#include "solar/relativity/units.h"

#include "solar/constants.h"

#include <cmath>
#include <stdexcept>

namespace solar::relativity {
namespace {

constexpr double solar_mass_kg = 1.98847e30;
constexpr double cubic_kilometres_to_cubic_metres = 1.0e9;
constexpr double kilometres_to_metres = 1.0e3;
constexpr double G_si =
    constants::G * cubic_kilometres_to_cubic_metres;
constexpr double c_si = constants::C_LIGHT * kilometres_to_metres;

double require_finite(double value, const char* quantity) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(quantity);
    }
    return value;
}

} // namespace

GeometricUnits GeometricUnits::from_mass_kg(double input_mass_kg) {
    if (!std::isfinite(input_mass_kg) || input_mass_kg <= 0.0) {
        throw std::invalid_argument("mass_kg must be finite and positive");
    }

    const double length_scale = G_si * input_mass_kg / (c_si * c_si);
    return {
        input_mass_kg,
        length_scale,
        length_scale / c_si,
    };
}

GeometricUnits GeometricUnits::from_solar_masses(double solar_masses) {
    if (!std::isfinite(solar_masses) || solar_masses <= 0.0) {
        throw std::invalid_argument(
            "solar_masses must be finite and positive");
    }
    return from_mass_kg(solar_masses * solar_mass_kg);
}

double GeometricUnits::length_si_to_M(double metres) const {
    return require_finite(metres, "metres must be finite") / M_length_m;
}

double GeometricUnits::length_M_to_si(double value_M) const {
    return require_finite(value_M, "length in M must be finite") * M_length_m;
}

double GeometricUnits::time_si_to_M(double seconds) const {
    return require_finite(seconds, "seconds must be finite") / M_time_s;
}

double GeometricUnits::time_M_to_si(double value_M) const {
    return require_finite(value_M, "time in M must be finite") * M_time_s;
}

double GeometricUnits::velocity_si_to_c(double metres_per_second) const {
    return require_finite(
               metres_per_second, "metres_per_second must be finite") /
           c_si;
}

} // namespace solar::relativity
