#include "solar/relativity/geodesic_types.h"
#include "solar/relativity/types.h"

#include <cmath>
#include <iostream>
#include <type_traits>

using namespace solar::relativity;

static_assert(!std::is_convertible_v<Contravariant4, Covariant4>);
static_assert(!std::is_convertible_v<Covariant4, Contravariant4>);

int main() {
    int passed = 0;
    int failed = 0;

    const auto check = [&passed, &failed](const char* name, bool condition) {
        std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
        condition ? ++passed : ++failed;
    };

    const PhaseSpaceState state{
        2.5,
        Contravariant4{Vec4{{1.0, 8.0, 0.5, -0.2}}},
        Covariant4{Vec4{{-0.9, 0.1, 2.0, 3.0}}},
    };
    check("affine parameter retained", state.affine == 2.5);
    check("coordinate order retained", state.x.v[0] == 1.0);
    check("covariant momentum retained", state.p.v[0] == -0.9);

    const GeodesicSample null_sample{};
    check("null proper time defaults unavailable",
          std::isnan(null_sample.proper_time));
    check("sample flags default clear", null_sample.flags == 0);
    check("null Hamiltonian target", hamiltonian_target(GeodesicKind::Null) == 0.0);
    check("timelike Hamiltonian target",
          hamiltonian_target(GeodesicKind::TimelikeUnitMass) == -0.5);
    const auto unknown_kind = static_cast<GeodesicKind>(99);
    check("unknown Hamiltonian target is unavailable",
          std::isnan(hamiltonian_target(unknown_kind)));
    check("proper time exists only for timelike",
          has_proper_time(GeodesicKind::TimelikeUnitMass) &&
          !has_proper_time(GeodesicKind::Null) &&
          !has_proper_time(unknown_kind));

    const IntegrationDiagnostics phase_one_aggregate{
        1,
        2,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        TerminationReason::MaxSteps,
        "phase-one aggregate",
    };
    check(
        "Phase 1 diagnostic aggregate reason remains compatible",
        phase_one_aggregate.reason ==
            TerminationReason::MaxSteps);
    check(
        "Phase 1 diagnostic aggregate message remains compatible",
        phase_one_aggregate.message ==
            "phase-one aggregate");
    check(
        "new absolute Carter diagnostic defaults unavailable",
        std::isnan(
            phase_one_aggregate.max_carter_abs_error));
    check(
        "near-critical termination reason is representable",
        TerminationReason::NearCriticalOrbit !=
            TerminationReason::UserEvent);

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
