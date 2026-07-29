#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace solar::relativity;

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

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
    double initial_step,
    double max_step,
    double max_affine) {
    return GeodesicIntegrationConfig::cpu_reference(
        GeodesicKind::Null,
        1.0,
        initial_step,
        max_step,
        max_affine);
}

} // namespace

int main() {
    const SchwarzschildBoyerLindquistMetric metric(1.0);
    const GeodesicIntegrator integrator(metric);

    const PhaseSpaceState outgoing{
        0.0,
        Contravariant4{Vec4{{0.0, 10.0, pi / 2.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.25, 0.0, 0.0}}},
    };
    const auto radial = integrator.integrate(
        outgoing, null_config(0.05, 0.1, 2.0));
    const double expected_time =
        2.0 + 2.0 * std::log((12.0 - 2.0) / (10.0 - 2.0));
    check("radial null reaches affine limit",
          radial.diagnostics.reason ==
              TerminationReason::MaxAffine);
    check_near("radial null radius",
               radial.final_state.x.v[1], 12.0, 2.0e-10);
    check_near("radial null coordinate time",
               radial.final_state.x.v[0],
               expected_time, 2.0e-10);
    check("radial null Hamiltonian gate",
          radial.diagnostics.max_constraint_error < 1.0e-10);

    const PhaseSpaceState photon_sphere{
        0.0,
        Contravariant4{Vec4{{0.0, 3.0, pi / 2.0, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0, 3.0 * std::sqrt(3.0)}}},
    };
    const auto sphere_rhs =
        HamiltonGeodesicRhs(metric)(photon_sphere);
    check_near("photon-sphere Hamiltonian",
               hamiltonian(metric, photon_sphere),
               0.0, 2.0e-14);
    check_near("photon-sphere radial velocity",
               sphere_rhs.dx.v[1], 0.0, 2.0e-14);
    check_near("photon-sphere radial momentum derivative",
               sphere_rhs.dp.v[1], 0.0, 2.0e-14);

    constexpr double impact_parameter = 100.0;
    constexpr double endpoint_radius = 10000.0;
    PhaseSpaceState weak_field{
        0.0,
        Contravariant4{
            Vec4{{0.0, endpoint_radius, pi / 2.0, 0.0}}},
        Covariant4{
            Vec4{{-1.0, 0.0, 0.0, impact_parameter}}},
    };
    const Mat4 inverse = metric.contravariant(weak_field.x);
    const double nonradial =
        inverse[0][0] -
        2.0 * impact_parameter * inverse[0][3] +
        impact_parameter * impact_parameter * inverse[3][3];
    weak_field.p.v[1] =
        -std::sqrt(-nonradial / inverse[1][1]);

    const GeodesicEvent escaped{
        "return to weak-field radius",
        [](const PhaseSpaceState& state) {
            return state.x.v[1] - endpoint_radius;
        },
        EventDirection::Increasing,
        TerminationReason::Escaped,
        1.0e-10,
    };
    const auto bending = integrator.integrate(
        weak_field,
        null_config(1.0, 10.0, 30000.0),
        {escaped});
    const double flat_finite_angle =
        pi - 2.0 * std::asin(
                 impact_parameter / endpoint_radius);
    const double deflection =
        std::fabs(
            bending.final_state.x.v[3] -
            weak_field.x.v[3]) -
        flat_finite_angle;
    const double leading_order =
        4.0 / impact_parameter;

    check("weak-field ray escapes",
          bending.diagnostics.reason ==
              TerminationReason::Escaped);
    check("weak-field event payload exists",
          bending.event.has_value());
    check_near("weak-field deflection approaches 4M/b",
               deflection, leading_order,
               0.05 * leading_order);
    check("weak-field Hamiltonian gate",
          bending.diagnostics.max_constraint_error < 1.0e-10);

    std::cout << std::setprecision(17)
              << "  radial_max_constraint="
              << radial.diagnostics.max_constraint_error
              << " weak_deflection=" << deflection
              << " weak_target=" << leading_order
              << " weak_max_constraint="
              << bending.diagnostics.max_constraint_error
              << '\n';
    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
