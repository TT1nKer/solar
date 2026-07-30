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
    check("proper time exists only for timelike",
          has_proper_time(GeodesicKind::TimelikeUnitMass) &&
          !has_proper_time(GeodesicKind::Null));

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
