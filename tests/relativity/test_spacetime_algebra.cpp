#include "solar/relativity/spacetime_algebra.h"

#include <cmath>
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
    const Mat4 minkowski{{
        {{-1.0, 0.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0, 0.0}},
        {{0.0, 0.0, 1.0, 0.0}},
        {{0.0, 0.0, 0.0, 1.0}},
    }};
    const Contravariant4 u{
        Vec4{{2.0, 0.5, -1.0, 3.0}}};
    const Contravariant4 v{
        Vec4{{1.5, -2.0, 4.0, 0.25}}};

    check_near(
        "literal Minkowski inner product",
        metric_inner_product(minkowski, u, v),
        -7.25,
        0.0);

    const Covariant4 lowered = lower_index(minkowski, u);
    check_near(
        "lowered time component",
        lowered.v[0],
        -2.0,
        0.0);
    check_near(
        "lowered spatial component",
        lowered.v[2],
        -1.0,
        0.0);

    const Contravariant4 raised =
        raise_index(minkowski, lowered);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near(
            "Minkowski raise/lower round trip",
            raised.v[component],
            u.v[component],
            0.0);
    }
    check_near(
        "covector-vector pairing",
        covector_vector_pairing(lowered, v),
        -7.25,
        0.0);

    const Mat4 mixed_metric{{
        {{-2.0, 0.5, 0.0, 0.0}},
        {{0.5, 3.0, 0.0, 0.0}},
        {{0.0, 0.0, 2.0, 0.0}},
        {{0.0, 0.0, 0.0, 4.0}},
    }};
    const Mat4 mixed_inverse{{
        {{-0.48, 0.08, 0.0, 0.0}},
        {{0.08, 0.32, 0.0, 0.0}},
        {{0.0, 0.0, 0.5, 0.0}},
        {{0.0, 0.0, 0.0, 0.25}},
    }};
    const Contravariant4 mixed_u{
        Vec4{{1.0, 2.0, 3.0, 4.0}}};
    const Contravariant4 mixed_v{
        Vec4{{-1.0, 0.5, -2.0, 1.0}}};
    check_near(
        "non-diagonal metric inner product",
        metric_inner_product(
            mixed_metric, mixed_u, mixed_v),
        8.25,
        0.0);
    const Covariant4 mixed_lowered =
        lower_index(mixed_metric, mixed_u);
    check_near(
        "non-diagonal lowering includes cross term",
        mixed_lowered.v[0],
        -1.0,
        0.0);
    const Contravariant4 mixed_raised =
        raise_index(mixed_inverse, mixed_lowered);
    for (std::size_t component = 0; component < 4; ++component) {
        check_near(
            "non-diagonal raise/lower round trip",
            mixed_raised.v[component],
            mixed_u.v[component],
            1.0e-15);
    }

    Contravariant4 non_finite = u;
    non_finite.v[1] =
        std::numeric_limits<double>::quiet_NaN();
    check(
        "non-finite contraction remains non-finite",
        std::isnan(metric_inner_product(
            minkowski, non_finite, v)));

    std::cout << "\n=== Results: " << passed
              << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
