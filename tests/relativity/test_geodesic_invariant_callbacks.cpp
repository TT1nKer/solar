#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/minkowski_metric.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ")
              << name << '\n';
    condition ? ++passed : ++failed;
}

void check_near(
    const char* name,
    double actual,
    double expected,
    double tolerance) {
    check(
        name,
        std::isfinite(actual) &&
            std::fabs(actual - expected) <= tolerance);
}

GeodesicIntegrationConfig config() {
    return GeodesicIntegrationConfig::cpu_reference(
        GeodesicKind::Null, 1.0, 0.1, 0.2, 1.0);
}

PhaseSpaceState photon() {
    return PhaseSpaceState{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 0.6, 0.8, 0.0}}},
    };
}

bool contains(
    const std::string& text,
    const std::string& fragment) {
    return text.find(fragment) != std::string::npos;
}

} // namespace

int main() {
    const MinkowskiMetric metric;
    const GeodesicIntegrator integrator(metric);

    auto custom_config = config();
    custom_config.monitor_energy = true;
    custom_config.monitor_lz = true;
    custom_config.stationary_energy_evaluator =
        [](const PhaseSpaceState& state) {
            return state.x.v[1] + state.p.v[0];
        };
    custom_config.axial_angular_momentum_evaluator =
        [](const PhaseSpaceState& state) {
            return state.x.v[2] - 2.0 * state.p.v[3];
        };
    const GeodesicIntegrationResult custom =
        integrator.integrate(photon(), custom_config);

    check("custom evaluator integration reaches limit",
          custom.diagnostics.reason ==
              TerminationReason::MaxAffine);
    check_near("custom stationary energy evaluator is used",
               custom.diagnostics.max_energy_rel_error,
               0.6, 2.0e-14);
    check_near("custom axial angular momentum evaluator is used",
               custom.diagnostics.max_lz_rel_error,
               0.8, 2.0e-14);

    auto default_config = config();
    default_config.monitor_energy = true;
    default_config.monitor_lz = true;
    const GeodesicIntegrationResult defaults =
        integrator.integrate(photon(), default_config);
    check_near("empty energy evaluator preserves BL default",
               defaults.diagnostics.max_energy_rel_error,
               0.0, 0.0);
    check_near("empty Lz evaluator preserves BL default",
               defaults.diagnostics.max_lz_rel_error,
               0.0, 0.0);

    int disabled_calls = 0;
    auto disabled_config = config();
    disabled_config.stationary_energy_evaluator =
        [&disabled_calls](const PhaseSpaceState&) {
            ++disabled_calls;
            return 0.0;
        };
    disabled_config.axial_angular_momentum_evaluator =
        [&disabled_calls](const PhaseSpaceState&) {
            ++disabled_calls;
            return 0.0;
        };
    const GeodesicIntegrationResult disabled =
        integrator.integrate(photon(), disabled_config);
    check("disabled evaluators do not run",
          disabled_calls == 0);
    check("disabled evaluator integration reaches limit",
          disabled.diagnostics.reason ==
              TerminationReason::MaxAffine);

    auto throwing_energy_config = config();
    throwing_energy_config.monitor_energy = true;
    throwing_energy_config.stationary_energy_evaluator =
        [](const PhaseSpaceState&) -> double {
            throw std::runtime_error(
                "synthetic evaluator failure");
        };
    const GeodesicIntegrationResult throwing_energy =
        integrator.integrate(photon(), throwing_energy_config);
    check("throwing energy evaluator fails explicitly",
          throwing_energy.diagnostics.reason ==
              TerminationReason::NonFiniteState);
    check("throwing energy failure identifies invariant",
          contains(
              throwing_energy.diagnostics.message,
              "stationary energy"));
    check("throwing energy evaluator accepts no step",
          throwing_energy.diagnostics.accepted_steps == 0);

    auto non_finite_energy_config = config();
    non_finite_energy_config.monitor_energy = true;
    non_finite_energy_config.stationary_energy_evaluator =
        [](const PhaseSpaceState&) {
            return std::numeric_limits<double>::quiet_NaN();
        };
    const GeodesicIntegrationResult non_finite_energy =
        integrator.integrate(photon(), non_finite_energy_config);
    check("non-finite energy evaluator fails explicitly",
          non_finite_energy.diagnostics.reason ==
              TerminationReason::NonFiniteState);
    check("non-finite energy failure identifies invariant",
          contains(
              non_finite_energy.diagnostics.message,
              "stationary energy"));

    auto throwing_lz_config = config();
    throwing_lz_config.monitor_lz = true;
    throwing_lz_config.axial_angular_momentum_evaluator =
        [](const PhaseSpaceState&) -> double {
            throw std::runtime_error(
                "synthetic evaluator failure");
        };
    const GeodesicIntegrationResult throwing_lz =
        integrator.integrate(photon(), throwing_lz_config);
    check("throwing Lz evaluator fails explicitly",
          throwing_lz.diagnostics.reason ==
              TerminationReason::NonFiniteState);
    check("throwing Lz failure identifies invariant",
          contains(
              throwing_lz.diagnostics.message,
              "axial angular momentum"));

    auto non_finite_lz_config = config();
    non_finite_lz_config.monitor_lz = true;
    non_finite_lz_config.axial_angular_momentum_evaluator =
        [](const PhaseSpaceState&) {
            return std::numeric_limits<double>::infinity();
        };
    const GeodesicIntegrationResult non_finite_lz =
        integrator.integrate(photon(), non_finite_lz_config);
    check("non-finite Lz evaluator fails explicitly",
          non_finite_lz.diagnostics.reason ==
              TerminationReason::NonFiniteState);
    check("non-finite Lz failure identifies invariant",
          contains(
              non_finite_lz.diagnostics.message,
              "axial angular momentum"));

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
