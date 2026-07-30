#pragma once

#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/observer.h"

#include <optional>
#include <string>

namespace solar::relativity {

enum class OrbitSense {
    Prograde,
    Retrograde,
};

constexpr bool is_valid_orbit_sense(
    OrbitSense sense) noexcept {
    return sense == OrbitSense::Prograde ||
           sense == OrbitSense::Retrograde;
}

enum class CircularOrbitStability {
    Stable,
    Unstable,
};

struct CircularTimelikeOrbit {
    double radius;
    double angular_velocity;
    double specific_energy;
    double specific_lz;
    CircularOrbitStability stability;
};

struct CircularOrbitResult {
    std::optional<CircularTimelikeOrbit> orbit;
    std::string message;

    explicit operator bool() const noexcept {
        return orbit.has_value();
    }
};

double kerr_isco_radius(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense);

double kerr_equatorial_photon_radius(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense);

double kerr_marginally_bound_radius(
    const KerrBoyerLindquistMetric& metric,
    OrbitSense sense);

CircularOrbitResult
evaluate_equatorial_circular_timelike_orbit(
    const KerrBoyerLindquistMetric& metric,
    double radius,
    OrbitSense sense);

ObserverResult make_equatorial_circular_observer(
    const KerrBoyerLindquistMetric& metric,
    const Contravariant4& x,
    OrbitSense sense);

} // namespace solar::relativity
