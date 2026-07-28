#include "solar/relativity/schwarzschild.h"
#include "solar/constants.h"
#include <cmath>
#include <stdexcept>

namespace solar::relativity {

SchwarzschildSpacetime::SchwarzschildSpacetime(double mass_length_km)
    : mass_length_km_(mass_length_km) {
    if (!std::isfinite(mass_length_km_) || mass_length_km_ <= 0.0)
        throw std::invalid_argument("Schwarzschild mass length must be finite and positive");
}

SchwarzschildSpacetime SchwarzschildSpacetime::from_mass_kg(double mass_kg) {
    if (!std::isfinite(mass_kg) || mass_kg <= 0.0)
        throw std::invalid_argument("Black-hole mass must be finite and positive");
    double c_squared = constants::C_LIGHT * constants::C_LIGHT;
    return SchwarzschildSpacetime(constants::G * mass_kg / c_squared);
}

static void require_exterior_position(
    const SchwarzschildSpacetime& spacetime, const FourVector& position) {
    if (!std::isfinite(position.r) || position.r <= spacetime.horizon_radius())
        throw std::domain_error("Schwarzschild coordinates require r > 2M");
    if (!std::isfinite(position.theta) || position.theta <= 0.0 ||
        position.theta >= constants::PI)
        throw std::domain_error("Schwarzschild polar angle must lie in (0, pi)");
}

double metric_norm(const SchwarzschildSpacetime& spacetime, const GeodesicState& state) {
    require_exterior_position(spacetime, state.position);
    double mass = spacetime.mass_length();
    double radius = state.position.r;
    double sin_theta = std::sin(state.position.theta);
    double lapse = 1.0 - 2.0 * mass / radius;
    const auto& u = state.tangent;
    return -lapse * u.t * u.t + u.r * u.r / lapse +
           radius * radius * (u.theta * u.theta + sin_theta * sin_theta * u.phi * u.phi);
}

ConservedQuantities conserved_quantities(
    const SchwarzschildSpacetime& spacetime, const GeodesicState& state) {
    require_exterior_position(spacetime, state.position);
    double radius = state.position.r;
    double lapse = 1.0 - 2.0 * spacetime.mass_length() / radius;
    double sin_theta = std::sin(state.position.theta);
    return {
        lapse * state.tangent.t,
        radius * radius * sin_theta * sin_theta * state.tangent.phi,
    };
}

std::vector<double> schwarzschild_geodesic_derivative(
    const SchwarzschildSpacetime& spacetime, const std::vector<double>& state) {
    if (state.size() != 8)
        throw std::invalid_argument("Schwarzschild geodesic state must contain 8 values");

    double mass = spacetime.mass_length();
    double radius = state[1];
    double theta = state[2];
    if (!std::isfinite(radius) || radius <= spacetime.horizon_radius())
        throw std::domain_error("Schwarzschild geodesic left the exterior coordinate chart");

    double ut = state[4];
    double ur = state[5];
    double utheta = state[6];
    double uphi = state[7];
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    if (std::fabs(sin_theta) < 1e-14)
        throw std::domain_error("Schwarzschild polar coordinate singularity reached");

    double radius_minus_2m = radius - 2.0 * mass;
    double radius_squared = radius * radius;
    double radius_cubed = radius_squared * radius;

    double dut = -2.0 * mass * ut * ur / (radius * radius_minus_2m);
    double dur = -mass * radius_minus_2m * ut * ut / radius_cubed
               + mass * ur * ur / (radius * radius_minus_2m)
               + radius_minus_2m * (utheta * utheta + sin_theta * sin_theta * uphi * uphi);
    double dutheta = -2.0 * ur * utheta / radius + sin_theta * cos_theta * uphi * uphi;
    double duphi = -2.0 * ur * uphi / radius
                 - 2.0 * cos_theta * utheta * uphi / sin_theta;

    return {ut, ur, utheta, uphi, dut, dur, dutheta, duphi};
}

GeodesicState make_equatorial_circular_timelike_state(
    const SchwarzschildSpacetime& spacetime, double radius_km) {
    double mass = spacetime.mass_length();
    if (!std::isfinite(radius_km) || radius_km <= 3.0 * mass)
        throw std::invalid_argument("Circular timelike orbit requires r > 3M");

    double normalization = std::sqrt(1.0 - 3.0 * mass / radius_km);
    GeodesicState state;
    state.position = {0.0, radius_km, constants::PI / 2.0, 0.0};
    state.tangent.t = 1.0 / normalization;
    state.tangent.phi = std::sqrt(mass / (radius_km * radius_km * radius_km)) / normalization;
    return state;
}

GeodesicState make_photon_sphere_state(const SchwarzschildSpacetime& spacetime) {
    double radius = 3.0 * spacetime.mass_length();
    double lapse = 1.0 - 2.0 * spacetime.mass_length() / radius;
    GeodesicState state;
    state.position = {0.0, radius, constants::PI / 2.0, 0.0};
    state.tangent.t = 1.0 / std::sqrt(lapse);
    state.tangent.phi = 1.0 / radius;
    return state;
}

} // namespace solar::relativity
