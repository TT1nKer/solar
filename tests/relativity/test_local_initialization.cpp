#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/local_initialization.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/observer.h"
#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace solar::relativity;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  PASS: " << name << "\n";
        ++passed;
    } else {
        std::cerr << "  FAIL: " << name << "\n";
        ++failed;
    }
}

void check_near(
    const std::string& name,
    double actual,
    double expected,
    double tolerance) {
    check(name, std::fabs(actual - expected) <= tolerance);
}

} // namespace

int main() {
    const MinkowskiMetric metric;
    const Contravariant4 position{
        Vec4{{0.0, 10.0, 0.0, 0.0}}};
    const ObserverResult observer =
        make_static_observer(metric, position);
    check("initialization observer exists", bool(observer));

    const InitialStateResult photon =
        initialize_local_photon(
            metric,
            *observer.frame,
            Vec3{{2.0, 0.0, 0.0}});
    check("local photon initializes", bool(photon));
    check(
        "successful photon has no error",
        photon.error == InitialStateError::None);
    check_near(
        "photon position retained",
        photon.state->x.v[1],
        10.0,
        0.0);
    check_near(
        "photon p_t",
        photon.state->p.v[0],
        -1.0,
        1.0e-15);
    check_near(
        "photon p_x",
        photon.state->p.v[1],
        1.0,
        1.0e-15);
    check_near(
        "photon observer frequency",
        photon.measured_frequency,
        1.0,
        1.0e-15);
    const double photon_constraint =
        hamiltonian_constraint_error(
            metric,
            *photon.state,
            GeodesicKind::Null);
    check(
        "photon null constraint",
        photon_constraint < 1.0e-14);

    const Contravariant4 photon_tangent = raise_index(
        metric.contravariant(position),
        photon.state->p);
    const Vec4 photon_local = coordinate_to_tetrad(
        metric.covariant(position),
        observer.frame->tetrad,
        photon_tangent);
    check_near(
        "photon local time round trip",
        photon_local[0],
        1.0,
        1.0e-15);
    check_near(
        "photon local direction round trip",
        photon_local[1],
        1.0,
        1.0e-15);
    check_near(
        "photon local transverse component",
        photon_local[2],
        0.0,
        0.0);

    const InitialStateResult timelike =
        initialize_local_timelike(
            metric,
            *observer.frame,
            Vec3{{0.6, 0.0, 0.0}});
    check("local timelike state initializes", bool(timelike));
    check_near(
        "timelike p_t includes gamma",
        timelike.state->p.v[0],
        -1.25,
        2.0e-15);
    check_near(
        "timelike p_x includes gamma v",
        timelike.state->p.v[1],
        0.75,
        2.0e-15);
    check_near(
        "timelike measured local energy",
        timelike.measured_frequency,
        1.25,
        2.0e-15);
    check_near(
        "timelike Hamiltonian target",
        hamiltonian(metric, *timelike.state),
        -0.5,
        2.0e-15);
    const double timelike_constraint =
        hamiltonian_constraint_error(
            metric,
            *timelike.state,
            GeodesicKind::TimelikeUnitMass);
    check(
        "timelike normalized constraint",
        timelike_constraint < 1.0e-14);

    check(
        "zero photon direction rejected",
        initialize_local_photon(
            metric,
            *observer.frame,
            Vec3{{0.0, 0.0, 0.0}}).error ==
            InitialStateError::InvalidLocalDirection);

    Vec3 non_finite_direction{{1.0, 0.0, 0.0}};
    non_finite_direction[1] =
        std::numeric_limits<double>::quiet_NaN();
    check(
        "non-finite photon direction rejected",
        initialize_local_photon(
            metric,
            *observer.frame,
            non_finite_direction).error ==
            InitialStateError::NonFiniteInput);

    check(
        "luminal local velocity rejected",
        initialize_local_timelike(
            metric,
            *observer.frame,
            Vec3{{1.0, 0.0, 0.0}}).error ==
            InitialStateError::SuperluminalLocalVelocity);
    check(
        "superluminal local velocity rejected",
        initialize_local_timelike(
            metric,
            *observer.frame,
            Vec3{{1.1, 0.0, 0.0}}).error ==
            InitialStateError::SuperluminalLocalVelocity);

    ObserverFrame invalid_observer = *observer.frame;
    invalid_observer.tetrad.basis[1].v =
        2.0 * invalid_observer.tetrad.basis[1].v;
    check(
        "non-orthonormal observer rejected",
        initialize_local_photon(
            metric,
            invalid_observer,
            Vec3{{1.0, 0.0, 0.0}}).error ==
            InitialStateError::InvalidObserverFrame);

    ObserverFrame outside_tetrad_gate = *observer.frame;
    outside_tetrad_gate.tetrad.basis[0].v[0] =
        std::sqrt(1.0 - 5.0e-11);
    check(
        "initialization rejects frame outside v3 tetrad gate",
        initialize_local_photon(
            metric,
            outside_tetrad_gate,
            Vec3{{1.0, 0.0, 0.0}}).error ==
            InitialStateError::InvalidObserverFrame);

    ObserverFrame inside_tetrad_gate = *observer.frame;
    inside_tetrad_gate.tetrad.basis[0].v[0] =
        std::sqrt(1.0 - 5.0e-13);
    check(
        "initialization accepts frame inside v3 tetrad gate",
        bool(initialize_local_photon(
            metric,
            inside_tetrad_gate,
            Vec3{{1.0, 0.0, 0.0}})));

    check(
        "non-finite affine parameter rejected",
        initialize_local_photon(
            metric,
            *observer.frame,
            Vec3{{1.0, 0.0, 0.0}},
            std::numeric_limits<double>::quiet_NaN()).error ==
            InitialStateError::NonFiniteInput);

    GeodesicIntegrationConfig backward_config =
        GeodesicIntegrationConfig::cpu_reference(
            GeodesicKind::Null,
            1.0,
            -0.1,
            0.1,
            1.0);
    const GeodesicIntegrationResult backward =
        GeodesicIntegrator(metric).integrate(
            *photon.state, backward_config);
    check(
        "backward photon reaches affine limit",
        backward.diagnostics.reason ==
            TerminationReason::MaxAffine);
    check_near(
        "negative affine step traces toward smaller x",
        backward.final_state.x.v[1],
        9.0,
        2.0e-14);
    check(
        "backward-traced momentum remains future directed",
        observer_measured_frequency(
            backward.final_state.p,
            observer.frame->tetrad.basis[0]) > 0.0);
    check_near(
        "backward-traced frequency remains normalized",
        observer_measured_frequency(
            backward.final_state.p,
            observer.frame->tetrad.basis[0]),
        1.0,
        1.0e-15);

    std::cout << std::setprecision(17)
              << "  photon_frequency="
              << photon.measured_frequency
              << " photon_constraint="
              << photon_constraint
              << " timelike_frequency="
              << timelike.measured_frequency
              << " timelike_constraint="
              << timelike_constraint
              << " backward_frequency="
              << observer_measured_frequency(
                     backward.final_state.p,
                     observer.frame->tetrad.basis[0])
              << "\n";
    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
