#pragma once

#include "solar/relativity/types.h"
#include <vector>

namespace solar::relativity {

class SchwarzschildSpacetime {
public:
    explicit SchwarzschildSpacetime(double mass_length_km);

    static SchwarzschildSpacetime from_mass_kg(double mass_kg);

    double mass_length() const { return mass_length_km_; }
    double horizon_radius() const { return 2.0 * mass_length_km_; }

private:
    double mass_length_km_;
};

double metric_norm(const SchwarzschildSpacetime& spacetime, const GeodesicState& state);

ConservedQuantities conserved_quantities(
    const SchwarzschildSpacetime& spacetime, const GeodesicState& state);

// Returns d(position,tangent)/d(lambda) in the flat order used by the ODE solver.
std::vector<double> schwarzschild_geodesic_derivative(
    const SchwarzschildSpacetime& spacetime, const std::vector<double>& state);

GeodesicState make_equatorial_circular_timelike_state(
    const SchwarzschildSpacetime& spacetime, double radius_km);

GeodesicState make_photon_sphere_state(const SchwarzschildSpacetime& spacetime);

} // namespace solar::relativity
