#include "solar/relativity/hamiltonian.h"
#include "solar/relativity/kerr_bl_metric.h"
#include "solar/relativity/minkowski_metric.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <cmath>
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

template <typename Action>
void check_domain_error(const char* name, Action action) {
    try {
        action();
        check(name, false);
    } catch (const std::domain_error&) {
        check(name, true);
    }
}

class TimeDependentMetric final : public Metric {
public:
    Chart chart() const noexcept override {
        return Chart::MinkowskiCartesian;
    }

    std::string name() const override {
        return "time-dependent-test";
    }

    Mat4 covariant(const Contravariant4& x) const override {
        require_valid(x);
        Mat4 metric{};
        metric[0][0] = -1.0 / (1.0 + 0.1 * x.v[0]);
        metric[1][1] = 1.0;
        metric[2][2] = 1.0;
        metric[3][3] = 1.0;
        return metric;
    }

    Mat4 contravariant(const Contravariant4& x) const override {
        require_valid(x);
        Mat4 inverse{};
        inverse[0][0] = -(1.0 + 0.1 * x.v[0]);
        inverse[1][1] = 1.0;
        inverse[2][2] = 1.0;
        inverse[3][3] = 1.0;
        return inverse;
    }

    std::array<Mat4, 4> contravariant_derivatives(
        const Contravariant4& x) const override {
        require_valid(x);
        std::array<Mat4, 4> derivatives{};
        derivatives[0][0][0] = -0.1;
        return derivatives;
    }

    bool valid_point(const Contravariant4& x) const noexcept override {
        return x.v.all_finite() &&
               std::isfinite(1.0 + 0.1 * x.v[0]) &&
               1.0 + 0.1 * x.v[0] != 0.0;
    }

private:
    void require_valid(const Contravariant4& x) const {
        if (!valid_point(x)) {
            throw std::domain_error(
                "invalid point for time-dependent test metric");
        }
    }
};

} // namespace

int main() {
    const MinkowskiMetric minkowski;
    const PhaseSpaceState photon{
        0.0,
        Contravariant4{Vec4{{0.0, 0.0, 0.0, 0.0}}},
        Covariant4{Vec4{{-1.0, 1.0, 0.0, 0.0}}},
    };
    check_near("null Hamiltonian",
               hamiltonian(minkowski, photon), 0.0, 0.0);
    const auto photon_rhs = HamiltonGeodesicRhs(minkowski)(photon);
    check_near("future time tangent",
               photon_rhs.dx.v[0], 1.0, 0.0);
    check_near("null spatial tangent",
               photon_rhs.dx.v[1], 1.0, 0.0);
    check_near("Minkowski momentum derivative zero",
               max_norm(photon_rhs.dp.v), 0.0, 0.0);

    const PhaseSpaceState massive{
        0.0,
        photon.x,
        Covariant4{Vec4{{-1.25, 0.75, 0.0, 0.0}}},
    };
    check_near("unit-mass timelike Hamiltonian",
               hamiltonian(minkowski, massive), -0.5, 0.0);
    check_near("timelike normalized constraint",
               hamiltonian_constraint_error(
                   minkowski, massive,
                   GeodesicKind::TimelikeUnitMass),
               0.0, 0.0);

    PhaseSpaceState perturbed = photon;
    perturbed.p.v = Vec4{{-2.0, 2.001, 0.0, 0.0}};
    check_near("hand-computed normalized null constraint",
               hamiltonian_constraint_error(
                   minkowski, perturbed, GeodesicKind::Null),
               0.00039993998401239666, 5.0e-17);

    const auto packed = pack_phase_space(perturbed);
    check_near("packed first coordinate",
               packed[0], perturbed.x.v[0], 0.0);
    check_near("packed final momentum",
               packed[7], perturbed.p.v[3], 0.0);
    const auto unpacked = unpack_phase_space(3.5, packed);
    check_near("unpacked affine", unpacked.affine, 3.5, 0.0);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near("coordinate pack round-trip",
                   unpacked.x.v[component],
                   perturbed.x.v[component], 0.0);
        check_near("momentum pack round-trip",
                   unpacked.p.v[component],
                   perturbed.p.v[component], 0.0);
    }

    const TimeDependentMetric time_dependent;
    PhaseSpaceState evolving = photon;
    evolving.p.v = Vec4{{-1.0, 0.0, 0.0, 0.0}};
    const auto evolving_rhs =
        HamiltonGeodesicRhs(time_dependent)(evolving);
    check_near("time-dependent metric changes p_t",
               evolving_rhs.dp.v[0], 0.05, 0.0);

    const SchwarzschildBoyerLindquistMetric schwarzschild(1.0);
    const KerrBoyerLindquistMetric kerr(1.0, 0.7);
    const PhaseSpaceState exterior{
        0.0,
        Contravariant4{Vec4{{0.3, 10.0, 1.2, -0.4}}},
        Covariant4{Vec4{{-1.0, 0.1, 0.2, 2.0}}},
    };
    const auto schwarzschild_rhs =
        HamiltonGeodesicRhs(schwarzschild)(exterior);
    check_near("Schwarzschild stationary p_t",
               schwarzschild_rhs.dp.v[0], 0.0, 0.0);
    check_near("Schwarzschild axisymmetric p_phi",
               schwarzschild_rhs.dp.v[3], 0.0, 0.0);
    const auto kerr_rhs = HamiltonGeodesicRhs(kerr)(exterior);
    check_near("Kerr stationary p_t",
               kerr_rhs.dp.v[0], 0.0, 0.0);
    check_near("Kerr axisymmetric p_phi",
               kerr_rhs.dp.v[3], 0.0, 0.0);

    PhaseSpaceState invalid = photon;
    invalid.affine = std::numeric_limits<double>::quiet_NaN();
    check_domain_error("non-finite affine rejected", [&] {
        (void)hamiltonian(minkowski, invalid);
    });
    invalid = photon;
    invalid.p.v[2] = std::numeric_limits<double>::infinity();
    check_domain_error("non-finite momentum rejected", [&] {
        (void)HamiltonGeodesicRhs(minkowski)(invalid);
    });
    PhaseSpaceState overflowing_constraint = photon;
    overflowing_constraint.p.v =
        Vec4{{-1.0e154, 1.0e154, 0.0, 0.0}};
    check_domain_error("overflowing constraint scale rejected", [&] {
        (void)hamiltonian_constraint_error(
            minkowski,
            overflowing_constraint,
            GeodesicKind::Null);
    });

    std::cout << "\n=== Results: " << passed << " passed, "
              << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
