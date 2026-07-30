#pragma once

namespace solar::relativity {

/**
 * Converts SI quantities to and from mass-scaled geometrized units.
 *
 * Equations:
 *   M_length = G M_kg / c^2
 *   M_time   = G M_kg / c^3
 *
 * Units:
 *   Inputs and stored scales are SI. Dimensionless outputs are measured in M
 *   with G=c=1.
 *
 * Failure:
 *   Construction rejects non-finite or non-positive mass. Conversion methods
 *   reject non-finite inputs.
 *
 * Validation:
 *   tests/relativity/test_units.cpp
 */
struct GeometricUnits {
    double mass_kg = 0.0;
    double M_length_m = 0.0;
    double M_time_s = 0.0;

    static GeometricUnits from_mass_kg(double mass_kg);
    static GeometricUnits from_solar_masses(double solar_masses);

    double length_si_to_M(double metres) const;
    double length_M_to_si(double value_M) const;
    double time_si_to_M(double seconds) const;
    double time_M_to_si(double value_M) const;
    double velocity_si_to_c(double metres_per_second) const;
};

} // namespace solar::relativity
