#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/minkowski_metric.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(const char* name, double actual, double expected,
                double tolerance) {
    check(name, std::isfinite(actual) &&
                    std::fabs(actual - expected) <= tolerance);
}

GeodesicIntegrationConfig null_config(
    double max_affine) {
    return GeodesicIntegrationConfig::cpu_reference(
        GeodesicKind::Null, 1.0, 0.05, 0.5, max_affine);
}

} // namespace

int main() {
    const KerrBoyerLindquistMetric kerr(1.0, 0.7);
    const GeodesicIntegrator kerr_integrator(kerr);

    constexpr double p_t = -1.0;
    constexpr double p_phi = 2.0;
    PhaseSpaceState ordinary{
        0.0,
        Contravariant4{Vec4{{0.0, 10.0, 1.2, 0.0}}},
        Covariant4{Vec4{{p_t, 0.0, 0.0, p_phi}}},
    };
    const Mat4 inverse = kerr.contravariant(ordinary.x);
    const double nonradial =
        inverse[0][0] * p_t * p_t +
        2.0 * inverse[0][3] * p_t * p_phi +
        inverse[3][3] * p_phi * p_phi;
    ordinary.p.v[1] =
        -std::sqrt(-nonradial / inverse[1][1]);

    auto monitored_config = null_config(2.0);
    monitored_config.monitor_energy = true;
    monitored_config.monitor_lz = true;
    const auto monitored = kerr_integrator.integrate(
        ordinary, monitored_config);

    check_near("ordinary Kerr initial Hamiltonian",
               hamiltonian(kerr, ordinary), 0.0, 2.0e-15);
    check("ordinary Kerr reaches affine limit",
          monitored.diagnostics.reason ==
              TerminationReason::MaxAffine);
    check_near("Kerr monitored energy drift",
               monitored.diagnostics.max_energy_rel_error,
               0.0, 0.0);
    check_near("Kerr monitored Lz drift",
               monitored.diagnostics.max_lz_rel_error,
               0.0, 0.0);
    check("Kerr Carter diagnostic remains unavailable",
          std::isnan(
              monitored.diagnostics.max_carter_rel_error));
    check("ordinary Kerr Hamiltonian gate",
          monitored.diagnostics.max_constraint_error < 1.0e-10);

    const MinkowskiMetric minkowski;
    const GeodesicIntegrator minkowski_integrator(minkowski);
    const PhaseSpaceState photon{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };

    const auto unmonitored = minkowski_integrator.integrate(
        photon, null_config(1.0));
    check("unmonitored energy is unavailable",
          std::isnan(
              unmonitored.diagnostics.max_energy_rel_error));
    check("unmonitored Lz is unavailable",
          std::isnan(
              unmonitored.diagnostics.max_lz_rel_error));
    check("unmonitored Carter is unavailable",
          std::isnan(
              unmonitored.diagnostics.max_carter_rel_error));

    auto zero_lz_config = null_config(1.0);
    zero_lz_config.monitor_lz = true;
    const auto zero_lz = minkowski_integrator.integrate(
        photon, zero_lz_config);
    check_near("zero Lz uses absolute drift",
               zero_lz.diagnostics.max_lz_rel_error,
               0.0, 0.0);
    check("energy remains unavailable when only Lz monitored",
          std::isnan(
              zero_lz.diagnostics.max_energy_rel_error));

    std::cout << std::setprecision(17)
              << "  kerr_initial_H=" << hamiltonian(kerr, ordinary)
              << " kerr_max_constraint="
              << monitored.diagnostics.max_constraint_error
              << " accepted="
              << monitored.diagnostics.accepted_steps
              << " rejected="
              << monitored.diagnostics.rejected_steps
              << '\n';
    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
