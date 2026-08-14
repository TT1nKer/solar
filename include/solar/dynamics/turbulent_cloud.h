#pragma once

#include "solar/body.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace solar {
namespace dynamics {

// Turbulent molecular-cloud initial conditions, realized as a Gaussian
// random field with a prescribed power spectrum (the Kritsuk-type
// supersonic-isothermal setup: log-density spectral index -2, log-normal
// density PDF) plus solid-body rotation.
//
// The generator owns the Fourier mode sum, so both the particle sample and
// the continuous field can be queried; the test suite verifies the
// realized statistics (structure function, dispersion, rms velocity,
// angular momentum) directly against the configured targets.
class TurbulentCloudGenerator {
public:
    struct Config {
        double radius_pc = 1.0;        // cloud radius
        double mass_solar = 1.0e4;     // total gas mass
        std::size_t particle_count = 16384;
        double density_spectral_index = -2.0;  // P(k) ~ k^n for ln rho
        double sigma_ln_rho = 2.0;     // log-density dispersion
        double velocity_spectral_index = -2.0; // E(k) ~ k^n
        double mach_number = 8.0;      // v_rms = mach * c_s
        double sound_speed_km_s = 0.19;        // ~10 K molecular gas
        double rotation_omega_per_myr = 0.6;   // solid-body rotation rate
        std::size_t mode_shells = 24;  // log-spaced |k| shells
        std::size_t modes_per_shell = 84;      // ~2000 modes total
        std::uint64_t seed = 7;
        bool multi_mass = false;       // mass proportional to local density
    };

    struct Realization {
        std::vector<Body> particles;   // solar units (kg, km, km/s)
        double density_rms_error = 0.0;    // sampled vs target sigma_ln_rho
        double velocity_rms_km_s = 0.0;    // measured turbulent rms
        Vec3 angular_momentum{};           // measured J
        double mean_density_kg_km3 = 0.0;
    };

    // Expose the continuous fields for statistical verification.
    struct Fields {
        double log_density_at(const Vec3& position_km) const;
        Vec3 velocity_at(const Vec3& position_km) const;
    private:
        friend class TurbulentCloudGenerator;
        const TurbulentCloudGenerator* owner_ = nullptr;
    };

    explicit TurbulentCloudGenerator(Config config);
    ~TurbulentCloudGenerator();

    TurbulentCloudGenerator(const TurbulentCloudGenerator&) = delete;
    TurbulentCloudGenerator& operator=(const TurbulentCloudGenerator&) = delete;

    Realization generate() const;
    Fields fields() const;

private:
    struct Impl;
    Impl* impl_;
};

} // namespace dynamics
} // namespace solar
