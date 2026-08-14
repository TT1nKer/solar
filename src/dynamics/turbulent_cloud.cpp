#include "solar/dynamics/turbulent_cloud.h"

#include "solar/constants.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace solar {
namespace dynamics {
namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double parsec_km = 3.085677581e13;
constexpr double solar_mass_kg = 1.98892e30;
constexpr double myr_s = 365.25 * 86400.0 * 1.0e6;

struct Mode {
    Vec3 k{};    // wave vector (1/km)
    double amplitude = 0.0;
    Vec3 phase{}; // random phases (for vector fields: direction handled separately)
};

// A three-dimensional Gaussian random field via a discrete Fourier mode
// sum. Modes are drawn on log-spaced |k| shells between k_min = 2 pi / L
// and k_max = shells * k_min (about two decades), with random directions.
struct ModeSum {
    std::vector<Mode> modes;
    double k_min = 0.0;
    double normalization = 1.0;

    // Scalar field: f(x) = sum_k A_k cos(k.x + phi_k).
    double scalar(const Vec3& x) const {
        double value = 0.0;
        for (const Mode& mode : modes) {
            const double phase =
                mode.k.dot(x) + mode.phase.x; // phase stored in phase.x
            value += mode.amplitude * std::cos(phase);
        }
        return value * normalization;
    }

    // Solenoidal vector field: v(x) = sum_k A_k (n1 cos + n2 sin),
    // where n1,n2 are two random orthonormal directions perpendicular to k,
    // so div v = 0 by construction.
    Vec3 solenoidal(const Vec3& x) const {
        Vec3 value{};
        for (const Mode& mode : modes) {
            const double phase = mode.k.dot(x);
            const Vec3& n1 = mode.phase;         // reused storage: n1
            const Vec3 n2 = mode.k.cross(n1).normalized();
            value += (n1 * std::cos(phase) + n2 * std::sin(phase)) *
                     mode.amplitude;
        }
        return value * normalization;
    }
};

} // namespace

struct TurbulentCloudGenerator::Impl {
    Config config;
    ModeSum density_modes;
    ModeSum velocity_modes;

    // Cached scales.
    double radius_km = 0.0;
    double total_mass_kg = 0.0;
    double mean_density = 0.0;   // kg / km^3
    double sound_speed_km_s = 0.0;

    Impl(Config cfg) : config(cfg) {
        radius_km = config.radius_pc * parsec_km;
        total_mass_kg = config.mass_solar * solar_mass_kg;
        const double volume = (4.0 / 3.0) * pi *
            radius_km * radius_km * radius_km;
        mean_density = total_mass_kg / volume;
        sound_speed_km_s = config.sound_speed_km_s;

        std::mt19937_64 generator(config.seed);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        const double k_min = 2.0 * pi / radius_km;
        const double k_max = k_min * 64.0;  // two decades of scales
        density_modes.k_min = k_min;
        velocity_modes.k_min = k_min;

        // Shell spacing: logarithmic.
        const std::size_t shells = config.mode_shells;
        const double log_ratio =
            std::log(k_max / k_min) / static_cast<double>(shells - 1);

        for (std::size_t shell = 0; shell < shells; ++shell) {
            const double k = k_min * std::exp(log_ratio * shell);
            // Modes sit on log-spaced shells with a constant count per
            // shell, so the sum approximates an integral over d(ln k):
            // sigma^2 ~ sum A_k^2 ~ int A(k)^2 d(ln k), while the power
            // spectrum P(k) ~ k^n gives sigma^2 = int 4 pi k^3 P(k)
            // d(ln k) ~ k^(n+3). Hence A(k) ~ k^((n+3)/2).
            const double density_amplitude = std::pow(
                k / k_min, (config.density_spectral_index + 3.0) / 2.0);
            const double velocity_amplitude = std::pow(
                k / k_min, (config.velocity_spectral_index + 3.0) / 2.0);
            for (std::size_t m = 0; m < config.modes_per_shell; ++m) {
                // Random direction on the unit sphere.
                const double cos_theta = 2.0 * unit(generator) - 1.0;
                const double phi = 2.0 * pi * unit(generator);
                const double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                Vec3 direction{
                    sin_theta * std::cos(phi),
                    sin_theta * std::sin(phi),
                    cos_theta};

                Mode dmode;
                dmode.k = direction * k;
                dmode.amplitude = density_amplitude;
                dmode.phase = {2.0 * pi * unit(generator), 0.0, 0.0};
                density_modes.modes.push_back(dmode);

                // Two random orthonormal directions perpendicular to k.
                Vec3 n1 = direction.cross(
                    std::fabs(direction.z) < 0.9
                        ? Vec3{0.0, 0.0, 1.0}
                        : Vec3{1.0, 0.0, 0.0})
                              .normalized();
                // Store the first solenoidal basis vector in the mode;
                // the second is recomputed inside solenoidal().
                Mode vmode;
                vmode.k = direction * k;
                vmode.amplitude = velocity_amplitude;
                vmode.phase = n1;
                velocity_modes.modes.push_back(vmode);
            }
        }

        // Normalize the density field to the requested log-density
        // dispersion and the velocity field to the requested rms, measured
        // on a deterministic probe set (same seed-derived stream).
        const std::size_t probe_count = 4096;
        const auto probe_point = [&](std::mt19937_64& engine) {
            const double r = radius_km * std::cbrt(unit(engine));
            const double ct = 2.0 * unit(engine) - 1.0;
            const double ph = 2.0 * pi * unit(engine);
            const double st = std::sqrt(1.0 - ct * ct);
            return Vec3{r * st * std::cos(ph), r * st * std::sin(ph), r * ct};
        };

        std::mt19937_64 density_probe(config.seed ^ 0x9e3779b97f4a7c15ULL);
        double variance_sum = 0.0;
        for (std::size_t i = 0; i < probe_count; ++i) {
            const Vec3 point = probe_point(density_probe);
            const double value = density_modes.scalar(point);
            variance_sum += value * value;
        }
        const double measured_sigma =
            std::sqrt(variance_sum / static_cast<double>(probe_count));
        if (measured_sigma > 1.0e-30) {
            density_modes.normalization = config.sigma_ln_rho / measured_sigma;
        }

        std::mt19937_64 velocity_probe(config.seed ^ 0x85ebca6bULL);
        variance_sum = 0.0;
        for (std::size_t i = 0; i < probe_count; ++i) {
            const Vec3 point = probe_point(velocity_probe);
            variance_sum += velocity_modes.solenoidal(point).norm_sq();
        }
        const double measured_vrms = std::sqrt(
            variance_sum / static_cast<double>(probe_count));
        if (measured_vrms > 1.0e-30) {
            velocity_modes.normalization =
                config.mach_number * config.sound_speed_km_s / measured_vrms;
        }
    }

    double log_density(const Vec3& position) const {
        return density_modes.scalar(position);
    }

    Vec3 turbulent_velocity(const Vec3& position) const {
        return velocity_modes.solenoidal(position);
    }
};

TurbulentCloudGenerator::TurbulentCloudGenerator(Config config)
    : impl_(new Impl(config)) {}

TurbulentCloudGenerator::~TurbulentCloudGenerator() { delete impl_; }

TurbulentCloudGenerator::Realization TurbulentCloudGenerator::generate() const {
    Impl& impl = *impl_;
    const Config& config = impl.config;
    Realization realization;

    std::mt19937_64 generator(config.seed ^ 0x243f6a8885a308d3ULL);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    // Rejection sampling of the density field inside the sphere.
    std::vector<Vec3> positions;
    positions.reserve(config.particle_count);
    std::vector<double> densities;
    densities.reserve(config.particle_count);

    // Find the field maximum for the acceptance bound.
    double density_max = -1.0e300;
    for (std::size_t probe = 0; probe < 8192; ++probe) {
        const double r = impl.radius_km * std::cbrt(unit(generator));
        const double ct = 2.0 * unit(generator) - 1.0;
        const double ph = 2.0 * pi * unit(generator);
        const double st = std::sqrt(1.0 - ct * ct);
        const Vec3 point{r * st * std::cos(ph), r * st * std::sin(ph), r * ct};
        density_max = std::max(density_max, impl.log_density(point));
    }

    while (positions.size() < config.particle_count) {
        const double r = impl.radius_km * std::cbrt(unit(generator));
        const double ct = 2.0 * unit(generator) - 1.0;
        const double ph = 2.0 * pi * unit(generator);
        const double st = std::sqrt(1.0 - ct * ct);
        const Vec3 point{r * st * std::cos(ph), r * st * std::sin(ph), r * ct};
        const double s = impl.log_density(point);
        const double accept = std::exp(s - density_max);
        if (unit(generator) >= accept) continue;
        positions.push_back(point);
        densities.push_back(impl.mean_density * std::exp(s));
    }

    // Masses: equal, or proportional to the local density.
    const double total_mass = impl.total_mass_kg;
    std::vector<double> masses(positions.size(), 0.0);
    if (config.multi_mass) {
        double density_sum = 0.0;
        for (const double density : densities) density_sum += density;
        for (std::size_t i = 0; i < masses.size(); ++i) {
            masses[i] = total_mass * densities[i] / density_sum;
        }
    } else {
        for (double& mass : masses) mass = total_mass / masses.size();
    }

    const double omega = config.rotation_omega_per_myr / myr_s;
    realization.particles.reserve(positions.size());
    double velocity_sum_sq = 0.0;
    Vec3 angular_momentum{};
    for (std::size_t i = 0; i < positions.size(); ++i) {
        Vec3 velocity = impl.turbulent_velocity(positions[i]);
        // Solid-body rotation around the z axis.
        velocity += Vec3{
            -omega * positions[i].y, omega * positions[i].x, 0.0};
        velocity_sum_sq += velocity.norm_sq();

        Body body;
        body.name = "cloud";
        body.mass = masses[i];
        body.mu = constants::G * masses[i];
        body.state.pos = positions[i];
        body.state.vel = velocity;
        angular_momentum +=
            positions[i].cross(velocity) * masses[i];
        realization.particles.push_back(body);
    }

    realization.velocity_rms_km_s =
        std::sqrt(velocity_sum_sq / positions.size());
    realization.angular_momentum = angular_momentum;
    realization.mean_density_kg_km3 = impl.mean_density;

    // Measured log-density dispersion of the realized particle sample.
    double log_sum = 0.0;
    double log_sq_sum = 0.0;
    for (const double density : densities) {
        const double s = std::log(density / impl.mean_density);
        log_sum += s;
        log_sq_sum += s * s;
    }
    const std::size_t n = densities.size();
    const double mean_log = log_sum / n;
    const double sigma_log = std::sqrt(
        std::max(0.0, log_sq_sum / n - mean_log * mean_log));
    realization.density_rms_error = sigma_log;

    return realization;
}

TurbulentCloudGenerator::Fields TurbulentCloudGenerator::fields() const {
    Fields fields;
    fields.owner_ = this;
    return fields;
}

double TurbulentCloudGenerator::Fields::log_density_at(
    const Vec3& position_km) const {
    return owner_->impl_->log_density(position_km);
}

Vec3 TurbulentCloudGenerator::Fields::velocity_at(
    const Vec3& position_km) const {
    return owner_->impl_->turbulent_velocity(position_km);
}

} // namespace dynamics
} // namespace solar
