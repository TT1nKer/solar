#include "solar/body.h"
#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_gravity.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using solar::Body;
using solar::Vec3;
using solar::dynamics::BarnesHutGravity;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

Body make_body(const Vec3& pos, double mass) {
    Body body;
    body.name = "p";
    body.mass = mass;
    body.mu = solar::constants::G * mass;
    body.state.pos = pos;
    body.state.vel = {};
    return body;
}

// Direct-sum reference with the same Plummer softening.
std::vector<Vec3> direct_acceleration(
    const std::vector<Body>& bodies, double softening_km) {
    std::vector<Vec3> acc(bodies.size(), Vec3{});
    const double softening_sq = softening_km * softening_km;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        for (std::size_t j = i + 1; j < bodies.size(); ++j) {
            Vec3 direction = bodies[j].state.pos - bodies[i].state.pos;
            const double inverse =
                1.0 / std::sqrt(direction.norm_sq() + softening_sq);
            const double scale = inverse * inverse * inverse;
            acc[i] += direction * (bodies[j].mu * scale);
            acc[j] -= direction * (bodies[i].mu * scale);
        }
    }
    return acc;
}

} // namespace

int main() {
    // --- Two-body exact check ---------------------------------------------
    {
        std::vector<Body> bodies;
        bodies.push_back(make_body({0.0, 0.0, 0.0}, 1.0e30));
        bodies.push_back(make_body({1.0e6, 0.0, 0.0}, 2.0e30));
        BarnesHutGravity gravity({0.7, 1000.0});
        std::vector<Vec3> acc(2, Vec3{});
        gravity.compute(bodies, 0.0, acc);
        // Plummer vector form: |a| = G M r / (r^2 + eps^2)^(3/2).
        const double r2 = 1.0e6 * 1.0e6 + 1000.0 * 1000.0;
        const double expected =
            solar::constants::G * 2.0e30 * 1.0e6 /
            (r2 * std::sqrt(r2));
        check("two-body force magnitude",
              std::abs(std::abs(acc[0].x) - expected) < 1.0e-12 * expected);
        check("two-body force direction",
              acc[0].x > 0.0 && acc[1].x < 0.0);
        check("two-body potential energy",
              std::abs(gravity.potential_energy(bodies, 0.0) +
                       solar::constants::G * 1.0e30 * 2.0e30 /
                           std::sqrt(1.0e12 + 1.0e6)) <
                  1.0e-6);
    }

    // --- Theta scan vs direct summation ------------------------------------
    // Clustered (cloud-like) distribution: four Gaussian clumps plus a dense
    // core, with a four-decade mass range. This is the regime the tree is
    // built for; a uniform cube would be the pathological worst case.
    {
        std::mt19937 generator(42);
        std::normal_distribution<double> gauss(0.0, 1.0);
        std::uniform_real_distribution<double> mass(1.0e28, 1.0e32);
        std::uniform_real_distribution<double> clump(0.0, 1.0);
        const double centers[5][3] = {
            {-8.0e6, -8.0e6, -8.0e6}, {8.0e6, -8.0e6, 8.0e6},
            {-8.0e6, 8.0e6, 8.0e6},   {8.0e6, 8.0e6, -8.0e6},
            {0.0, 0.0, 0.0}};
        const double widths[5] = {1.5e6, 1.5e6, 1.5e6, 1.5e6, 4.0e5};
        std::vector<Body> bodies;
        for (int i = 0; i < 600; ++i) {
            int c = static_cast<int>(clump(generator) * 5.0) % 5;
            bodies.push_back(make_body(
                {centers[c][0] + widths[c] * gauss(generator),
                 centers[c][1] + widths[c] * gauss(generator),
                 centers[c][2] + widths[c] * gauss(generator)},
                mass(generator)));
        }
        const std::vector<Vec3> reference = direct_acceleration(bodies, 100.0);

        double mean_field = 0.0;
        for (const Vec3& a : reference) mean_field += a.norm();
        mean_field /= static_cast<double>(reference.size());
        const double floor = 1.0e-8 * mean_field;

        std::vector<double> errors;
        for (const double theta : {0.15, 0.3, 0.5, 0.7}) {
            BarnesHutGravity gravity({theta, 100.0});
            std::vector<Vec3> acc(bodies.size(), Vec3{});
            gravity.compute(bodies, 0.0, acc);
            double max_relative = 0.0;
            for (std::size_t i = 0; i < bodies.size(); ++i) {
                const double magnitude =
                    (acc[i] - reference[i]).norm() /
                    std::max(reference[i].norm(), floor);
                max_relative = std::max(max_relative, magnitude);
            }
            errors.push_back(max_relative);
            std::cout << "  theta=" << theta
                      << " max_relative_acc_error=" << max_relative << '\n';
        }
        // Barnes-Hut truncation model: worst-case error ~ (s/d)^2 <=
        // theta^2, and the realized error must vanish as theta -> 0.
        // Monotonicity is checked across the sorted theta ladder (a node
        // flips from "descend" to "monopole" discretely, so max error can
        // jump between adjacent rungs; it must never decrease overall).
        check("accuracy non-decreasing across the theta ladder",
              errors[0] <= errors[1] && errors[1] <= errors[2] &&
                  errors[2] <= errors[3]);
        check("theta=0.15 error < 5e-3 (near-direct)", errors[0] < 5.0e-3);
        check("theta=0.3 error < 1e-2", errors[1] < 1.0e-2);
        check("theta=0.5 error < 1e-1", errors[2] < 1.0e-1);
        check("theta=0.7 error < 2.5e-1", errors[3] < 2.5e-1);

        BarnesHutGravity gravity({0.7, 100.0});
        double tree_pe = gravity.potential_energy(bodies, 0.0);
        double direct_pe = 0.0;
        const double softening_sq = 100.0 * 100.0;
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            for (std::size_t j = i + 1; j < bodies.size(); ++j) {
                direct_pe -= bodies[i].mu * bodies[j].mass /
                             std::sqrt((bodies[j].state.pos -
                                        bodies[i].state.pos)
                                           .norm_sq() +
                                       softening_sq);
            }
        }
        check("potential energy matches direct sum",
              std::abs(tree_pe - direct_pe) < 1.0e-9 * std::abs(direct_pe));
    }

    std::cout << (failures == 0 ? "PASS: barnes-hut gravity checks"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}
