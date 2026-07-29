#include "solar/relativity/math.h"

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
    check(name, std::fabs(actual - expected) <= tolerance);
}

} // namespace

int main() {
    const Vec4 a{{1.0, -2.0, 3.0, -4.0}};
    const Vec4 b{{0.5, 2.0, -1.0, 8.0}};

    const Vec4 sum = a + b;
    check_near("vector addition component", sum[3], 4.0, 0.0);
    check_near("vector subtraction component", (a - b)[1], -4.0, 0.0);
    check_near("left scalar multiplication", (2.0 * a)[2], 6.0, 0.0);
    check_near("right scalar multiplication", (a * 0.5)[3], -2.0, 0.0);
    check_near("scalar division", (a / 2.0)[1], -1.0, 0.0);
    check_near("max norm", max_norm(a), 4.0, 0.0);
    check("finite vector detected", a.all_finite());

    Vec4 non_finite = a;
    non_finite[2] = std::numeric_limits<double>::infinity();
    check("non-finite vector detected", !non_finite.all_finite());

    check_near("explicit Minkowski contraction",
               minkowski_dot_minus_plus_plus_plus(a, b), -39.5, 0.0);

    const Mat4 matrix{{
        {{4.0, 7.0, 2.0, 3.0}},
        {{0.0, 5.0, 0.0, 1.0}},
        {{0.0, 0.0, 3.0, 0.0}},
        {{0.0, 0.0, 0.0, 2.0}},
    }};
    const Mat4 matrix_inverse = inverse(matrix);
    const Mat4 product = multiply(matrix, matrix_inverse);
    double identity_error = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            identity_error = std::max(
                identity_error, std::fabs(product[row][column] - expected));
        }
    }
    check("matrix inverse recovers identity", identity_error < 1.0e-14);

    const Vec4 matrix_vector = multiply(matrix, Vec4{{1.0, 2.0, 3.0, 4.0}});
    check_near("matrix-vector multiplication", matrix_vector[0], 36.0, 0.0);
    check("finite matrix detected", all_finite(matrix));

    Mat4 singular = matrix;
    singular[2] = singular[1];
    try {
        (void)inverse(singular);
        check("singular matrix rejected", false);
    } catch (const std::domain_error&) {
        check("singular matrix rejected", true);
    }

    try {
        (void)(a / 0.0);
        check("zero vector divisor rejected", false);
    } catch (const std::domain_error&) {
        check("zero vector divisor rejected", true);
    }

    std::cout << "\n=== Results: " << passed << " passed, " << failed
              << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
