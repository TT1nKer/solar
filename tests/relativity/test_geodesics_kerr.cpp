#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/kerr_constants.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/observer.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

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
    check("Kerr Carter absolute diagnostic remains unavailable",
          std::isnan(
              monitored.diagnostics.max_carter_abs_error));
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
    check("unmonitored Carter absolute error is unavailable",
          std::isnan(
              unmonitored.diagnostics.max_carter_abs_error));

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

    const Contravariant4 generic_position{
        Vec4{{0.0, 8.0, 1.1, 0.0}}};
    const ObserverResult generic_zamo =
        make_zamo_observer(kerr, generic_position);
    check("generic Kerr ray ZAMO exists", bool(generic_zamo));
    const InitialStateResult generic_photon =
        initialize_local_photon(
            kerr,
            *generic_zamo.frame,
            Vec3{{0.3, 0.4, 0.5}});
    check("generic Kerr photon initializes", bool(generic_photon));

    auto carter_config =
        GeodesicIntegrationConfig::cpu_reference(
            GeodesicKind::Null,
            1.0,
            0.02,
            0.1,
            5.0);
    carter_config.monitor_energy = true;
    carter_config.monitor_lz = true;
    carter_config.carter_evaluator =
        [&kerr](const PhaseSpaceState& state) {
            return evaluate_kerr_constants(
                kerr, state, GeodesicKind::Null).Q;
        };
    const GeodesicIntegrationResult carter_monitored =
        kerr_integrator.integrate(
            *generic_photon.state, carter_config);
    check(
        "generic Kerr ray reaches affine limit",
        carter_monitored.diagnostics.reason ==
            TerminationReason::MaxAffine);
    check(
        "generic Kerr Hamiltonian gate",
        carter_monitored.diagnostics.max_constraint_error <
            1.0e-10);
    check(
        "generic Kerr Carter relative gate",
        std::isfinite(
            carter_monitored.diagnostics.max_carter_rel_error) &&
        carter_monitored.diagnostics.max_carter_rel_error <
            1.0e-10);
    check(
        "generic Kerr Carter absolute gate",
        std::isfinite(
            carter_monitored.diagnostics.max_carter_abs_error) &&
        carter_monitored.diagnostics.max_carter_abs_error <
            1.0e-10);

    auto synthetic_config = null_config(2.0);
    synthetic_config.carter_evaluator =
        [](const PhaseSpaceState& state) {
            return 1.0e-6 + 1.0e-9 * state.affine;
        };
    const GeodesicIntegrationResult synthetic =
        minkowski_integrator.integrate(
            photon, synthetic_config);
    check_near(
        "small invariant uses v3 max-one denominator",
        synthetic.diagnostics.max_carter_rel_error,
        2.0e-9,
        2.0e-17);
    check_near(
        "small invariant absolute error reported",
        synthetic.diagnostics.max_carter_abs_error,
        2.0e-9,
        2.0e-17);

    auto non_finite_monitor_config = null_config(1.0);
    non_finite_monitor_config.carter_evaluator =
        [](const PhaseSpaceState&) {
            return std::numeric_limits<double>::quiet_NaN();
        };
    const GeodesicIntegrationResult non_finite_monitor =
        minkowski_integrator.integrate(
            photon, non_finite_monitor_config);
    check(
        "non-finite invariant evaluator fails explicitly",
        non_finite_monitor.diagnostics.reason ==
            TerminationReason::NonFiniteState);
    check(
        "non-finite invariant evaluator accepts no step",
        non_finite_monitor.diagnostics.accepted_steps == 0);

    auto throwing_monitor_config = null_config(1.0);
    throwing_monitor_config.carter_evaluator =
        [](const PhaseSpaceState&) -> double {
            throw std::domain_error("test invariant failure");
        };
    const GeodesicIntegrationResult throwing_monitor =
        minkowski_integrator.integrate(
            photon, throwing_monitor_config);
    check(
        "throwing invariant evaluator fails explicitly",
        throwing_monitor.diagnostics.reason ==
            TerminationReason::NonFiniteState);
    check(
        "throwing invariant evaluator accepts no step",
        throwing_monitor.diagnostics.accepted_steps == 0);

    std::cout << std::setprecision(17)
              << "  kerr_initial_H=" << hamiltonian(kerr, ordinary)
              << " kerr_max_constraint="
              << monitored.diagnostics.max_constraint_error
              << " carter_rel="
              << carter_monitored.diagnostics.max_carter_rel_error
              << " carter_abs="
              << carter_monitored.diagnostics.max_carter_abs_error
              << " accepted="
              << monitored.diagnostics.accepted_steps
              << " rejected="
              << monitored.diagnostics.rejected_steps
              << '\n';
    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
