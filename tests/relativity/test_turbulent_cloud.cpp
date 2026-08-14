#include "solar/constants.h"
#include "solar/dynamics/turbulent_cloud.h"

#include <cmath>
#include <iostream>
#include <vector>

using solar::Vec3;
using solar::dynamics::TurbulentCloudGenerator;

namespace {


int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}


} // namespace

int main() {
    TurbulentCloudGenerator::Config config;
    config.radius_pc = 1.0;
    config.mass_solar = 1.0e4;
    config.particle_count = 16384;
    config.sigma_ln_rho = 2.0;
    config.mach_number = 8.0;
    config.sound_speed_km_s = 0.19;
    config.rotation_omega_per_myr = 0.6;
    config.seed = 7;

    const TurbulentCloudGenerator generator(config);
    const auto cloud = generator.generate();

    // 1. Determinism: same seed -> identical particles.
    const TurbulentCloudGenerator generator2(config);
    const auto cloud2 = generator2.generate();
    bool identical = cloud.particles.size() == cloud2.particles.size();
    for (std::size_t i = 0; identical && i < cloud.particles.size(); ++i) {
        identical = cloud.particles[i].state.pos.x ==
                        cloud2.particles[i].state.pos.x &&
                    cloud.particles[i].state.pos.y ==
                        cloud2.particles[i].state.pos.y &&
                    cloud.particles[i].state.pos.z ==
                        cloud2.particles[i].state.pos.z &&
                    cloud.particles[i].state.vel.x ==
                        cloud2.particles[i].state.vel.x &&
                    cloud.particles[i].mass == cloud2.particles[i].mass;
    }
    check("seeded generation is deterministic", identical);

    // 2. Total mass and particle count.
    double total_mass = 0.0;
    for (const auto& body : cloud.particles) total_mass += body.mass;
    check("total mass matches configuration",
          std::abs(total_mass / cloud.particles.size() /
                       (1.0e4 * 1.98892e30 / 16384.0) -
                   1.0) < 1.0e-12);
    check("particle count matches configuration",
          cloud.particles.size() == 16384);

    // 3. Realized log-density dispersion (particle sample) near target.
    //    Rejection sampling biases slightly high; allow 15% tolerance.
    std::cout << "  measured sigma_ln_rho = " << cloud.density_rms_error
              << " (target " << config.sigma_ln_rho << ")" << '\n';
    check("realized sigma_ln_rho within 15% of target",
          std::abs(cloud.density_rms_error - config.sigma_ln_rho) <
              0.15 * config.sigma_ln_rho);

    // 4. Velocity rms near mach * c_s (measured with rotation included,
    //    so allow a band).
    const double target_vrms = config.mach_number * config.sound_speed_km_s;
    const double measured = cloud.velocity_rms_km_s;
    std::cout << "  measured v_rms = " << measured
              << " km/s (target " << target_vrms << ")" << '\n';
    check("velocity rms within 25% of mach*c_s",
          std::abs(measured - target_vrms) < 0.25 * target_vrms);

    // 5. Angular momentum: the generator adds omega x r exactly, so
    //    J_z = omega * sum m_i (x_i^2 + y_i^2) must hold to machine
    //    precision for the realized (density-weighted) sample. The
    //    uniform-sphere value 2/5 M R^2 omega is reported as context
    //    only: the density field moves mass around, so it is not a gate.
    const double omega = config.rotation_omega_per_myr /
                         (365.25 * 86400.0 * 1.0e6);
    const double radius_km = config.radius_pc * parsec_km;
    double sampled_moment = 0.0;
    for (const auto& body : cloud.particles) {
        sampled_moment += body.mass *
            (body.state.pos.x * body.state.pos.x +
             body.state.pos.y * body.state.pos.y);
    }
    const double expected_jz = omega * sampled_moment;
    const double uniform_jz = 0.4 * (1.0e4 * 1.98892e30) *
                              radius_km * radius_km * omega;
    std::cout << "  J_z measured = " << cloud.angular_momentum.z
              << ", omega * I_sampled = " << expected_jz
              << ", uniform-sphere value = " << uniform_jz << '\n';
    // The turbulent velocity field also carries a small net angular
    // momentum (random walk over ~2k modes), so gate the rotation
    // identity at the 5% level rather than machine precision.
    check("rotation J_z equals omega times the sampled moment of inertia",
          std::abs(cloud.angular_momentum.z - expected_jz) <
              0.05 * std::abs(expected_jz));

    std::cout << (failures == 0 ? "PASS: turbulent cloud statistics"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}
