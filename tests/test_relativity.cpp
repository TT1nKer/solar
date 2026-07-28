#include "solar/relativity/geodesic.h"
#include "solar/constants.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace solar::relativity;

static int passed = 0;
static int failed = 0;

static void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << "\n";
    condition ? passed++ : failed++;
}

static GeodesicOptions accurate_options(double duration, double step) {
    GeodesicOptions options;
    options.affine_duration = duration;
    options.initial_step = step;
    options.absolute_tolerance = 1e-12;
    options.relative_tolerance = 1e-12;
    return options;
}

int main() {
    std::cout << "=== Schwarzschild Geodesic Tests ===\n";
    SchwarzschildSpacetime spacetime(1.0);

    {
        auto state = make_equatorial_circular_timelike_state(spacetime, 10.0);
        check("circular timelike initial state is normalized",
              std::fabs(metric_norm(spacetime, state) + 1.0) < 1e-14);

        double period = 2.0 * solar::constants::PI / state.tangent.phi;
        auto result = propagate_geodesic(
            spacetime, state, GeodesicKind::Timelike, accurate_options(period, 0.25));
        const auto& final_state = result.samples.back().state;
        check("r=10M circular orbit completes", result.status == PropagationStatus::Completed);
        check("r=10M radius remains circular",
              std::fabs(final_state.position.r - 10.0) < 1e-8);
        check("timelike metric norm is conserved", result.max_metric_norm_drift < 1e-10);
        check("timelike Killing invariants are conserved",
              result.max_relative_energy_drift < 1e-10 &&
              result.max_relative_angular_momentum_drift < 1e-10);
    }

    {
        auto state = make_photon_sphere_state(spacetime);
        check("photon-sphere initial tangent is null",
              std::fabs(metric_norm(spacetime, state)) < 1e-14);

        double period = 2.0 * solar::constants::PI / state.tangent.phi;
        auto result = propagate_geodesic(
            spacetime, state, GeodesicKind::Null, accurate_options(period, 0.05));
        const auto& final_state = result.samples.back().state;
        check("r=3M photon orbit completes", result.status == PropagationStatus::Completed);
        check("r=3M photon orbit remains circular",
              std::fabs(final_state.position.r - 3.0) < 1e-8);
        check("null metric norm is conserved", result.max_metric_norm_drift < 1e-10);
    }

    {
        auto solar_mass = SchwarzschildSpacetime::from_mass_kg(1.98847e30);
        check("solar mass converts to geometric length",
              solar_mass.mass_length() > 1.476 && solar_mass.mass_length() < 1.477);
    }

    {
        GeodesicState falling;
        falling.position = {0.0, 10.0, solar::constants::PI / 2.0, 0.0};
        falling.tangent.t = 1.0 / std::sqrt(1.0 - 2.0 / falling.position.r);
        auto result = propagate_geodesic(
            spacetime, falling, GeodesicKind::Timelike, accurate_options(5.0, 0.02));
        check("released particle falls inward",
              result.status == PropagationStatus::Completed &&
              result.samples.back().state.position.r < falling.position.r);
        check("infall conserves normalization and energy",
              result.max_metric_norm_drift < 1e-10 &&
              result.max_relative_energy_drift < 1e-10);

        auto capture_options = accurate_options(100.0, 0.02);
        capture_options.horizon_margin = 1e-4;
        auto capture = propagate_geodesic(
            spacetime, falling, GeodesicKind::Timelike, capture_options);
        check("radial infall terminates at the exterior horizon guard",
              capture.status == PropagationStatus::HorizonReached);
    }

    {
        bool rejected_horizon_state = false;
        try {
            GeodesicState invalid;
            invalid.position = {0.0, 2.0, solar::constants::PI / 2.0, 0.0};
            invalid.tangent.t = 1.0;
            propagate_geodesic(
                spacetime, invalid, GeodesicKind::Timelike, accurate_options(1.0, 0.1));
        } catch (const std::domain_error&) {
            rejected_horizon_state = true;
        }
        check("state at the Schwarzschild horizon is rejected", rejected_horizon_state);
    }

    std::cout << "=== Results: " << passed << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
