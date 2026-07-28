#include "blackhole_command.h"
#include "solar/constants.h"
#include "solar/relativity/geodesic.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

constexpr double SOLAR_MASS_KG = 1.98847e30;

double parse_positive_number(const char* text, const char* name) {
    size_t consumed = 0;
    double value = std::stod(text, &consumed);
    if (text[consumed] != '\0' || !std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument(std::string(name) + " must be a positive number");
    return value;
}

const char* status_name(solar::relativity::PropagationStatus status) {
    using solar::relativity::PropagationStatus;
    switch (status) {
        case PropagationStatus::Completed: return "completed";
        case PropagationStatus::HorizonReached: return "horizon_reached";
        case PropagationStatus::StepLimitReached: return "step_limit_reached";
        case PropagationStatus::StepSizeUnderflow: return "step_size_underflow";
    }
    return "unknown";
}

} // namespace

int run_blackhole_command(int argc, char* argv[]) {
    using namespace solar::relativity;

    if (argc < 5) {
        std::cerr << "Usage:\n"
                  << "  solar blackhole circular <mass_solar> <radius_M> <orbits>\n"
                  << "  solar blackhole photon   <mass_solar> <orbits>\n";
        return 1;
    }

    std::string mode = argv[2];
    double mass_solar = parse_positive_number(argv[3], "mass_solar");
    SchwarzschildSpacetime spacetime =
        SchwarzschildSpacetime::from_mass_kg(mass_solar * SOLAR_MASS_KG);

    GeodesicState initial_state;
    double orbit_count = 0.0;
    GeodesicKind kind;

    if (mode == "circular") {
        if (argc < 6) throw std::invalid_argument("circular mode requires radius_M and orbits");
        double radius_in_mass_lengths = parse_positive_number(argv[4], "radius_M");
        orbit_count = parse_positive_number(argv[5], "orbits");
        initial_state = make_equatorial_circular_timelike_state(
            spacetime, radius_in_mass_lengths * spacetime.mass_length());
        kind = GeodesicKind::Timelike;
    } else if (mode == "photon") {
        orbit_count = parse_positive_number(argv[4], "orbits");
        initial_state = make_photon_sphere_state(spacetime);
        kind = GeodesicKind::Null;
    } else {
        throw std::invalid_argument("blackhole mode must be 'circular' or 'photon'");
    }

    double affine_period = 2.0 * solar::constants::PI / initial_state.tangent.phi;
    GeodesicOptions options;
    options.affine_duration = orbit_count * affine_period;
    options.initial_step = affine_period / 1000.0;
    options.absolute_tolerance = 1e-11;
    options.relative_tolerance = 1e-11;

    GeodesicResult result = propagate_geodesic(spacetime, initial_state, kind, options);
    const GeodesicSample& final_sample = result.samples.back();

    std::cout << std::scientific << std::setprecision(10)
              << "Schwarzschild geodesic\n"
              << "  status: " << status_name(result.status) << "\n"
              << "  M: " << spacetime.mass_length() << " km\n"
              << "  horizon: " << spacetime.horizon_radius() << " km\n"
              << "  final r/M: "
              << final_sample.state.position.r / spacetime.mass_length() << "\n"
              << "  final phi: " << final_sample.state.position.phi << " rad\n"
              << "  accepted/rejected steps: " << result.accepted_steps
              << "/" << result.rejected_steps << "\n"
              << "  max metric-norm drift: " << result.max_metric_norm_drift << "\n"
              << "  max relative energy drift: " << result.max_relative_energy_drift << "\n"
              << "  max relative Lz drift: "
              << result.max_relative_angular_momentum_drift << "\n";

    return result.status == PropagationStatus::Completed ? 0 : 2;
}
