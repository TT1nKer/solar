#include "solar/body.h"
#include "solar/constants.h"
#include "solar/dynamics/barnes_hut_gravity.h"
#include "solar/dynamics/pn_gravity.h"
#include "solar/gravity.h"
#include "solar/relativity/geodesic_integrator.h"
#include "solar/relativity/schwarzschild_metric.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace solar;
using namespace solar::dynamics;
namespace gr = solar::relativity;

namespace {

int failures = 0;

void check(const std::string& name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << '\n';
    if (!condition) ++failures;
}

constexpr double pi = 3.14159265358979323846;

Body make_body(const Vec3& pos, const Vec3& vel, double mass) {
    Body body;
    body.name = "p";
    body.mass = mass;
    body.mu = constants::G * mass;
    body.state.pos = pos;
    body.state.vel = vel;
    return body;
}

} // namespace

int main() {
    const double c = constants::C_LIGHT;          // km/s
    const double mass_sun = 1.98892e30;           // kg
    const double mu_sun = constants::G * mass_sun;  // km^3/s^2
    const double rg = mu_sun / (c * c);           // GM/c^2 (km)

    // --- 1. c -> infinity reduces to Newtonian gravity ---------------------
    {
        std::mt19937 generator(11);
        std::uniform_real_distribution<double> position(-1.0e6, 1.0e6);
        std::uniform_real_distribution<double> velocity(-50.0, 50.0);
        std::uniform_real_distribution<double> mass(1.0e24, 1.0e28);
        std::vector<Body> bodies;
        for (int i = 0; i < 100; ++i) {
            bodies.push_back(make_body(
                {position(generator), position(generator),
                 position(generator)},
                {velocity(generator), velocity(generator),
                 velocity(generator)},
                mass(generator)));
        }
        NewtonianGravity newtonian;
        PostNewtonianGravity pn(c * 1.0e6);  // artificially huge c
        std::vector<Vec3> acc_newton(bodies.size(), Vec3{});
        std::vector<Vec3> acc_pn(bodies.size(), Vec3{});
        newtonian.compute(bodies, 0.0, acc_newton);
        pn.compute(bodies, 0.0, acc_pn);
        double max_relative = 0.0;
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            max_relative = std::max(
                max_relative,
                (acc_pn[i] - acc_newton[i]).norm() /
                    std::max(acc_newton[i].norm(), 1.0e-12));
        }
        std::cout << "  c->infinity max relative deviation = "
                  << max_relative << '\n';
        check("1PN force reduces to Newtonian as c -> infinity",
              max_relative < 1.0e-6);
    }

    // --- 2. Analytic radial-acceleration spot check -------------------------
    // At a reference point the 1PN radial acceleration must match the
    // standard gauge expression
    //   a_r = -(GM/r^2) [1 - 4 GM/(r c^2) - 3 v^2/c^2].
    {
        const double r0 = 200.0 * rg;
        Body central = make_body({0, 0, 0}, {0, 0, 0}, mass_sun);
        Body particle = make_body({r0, 0, 0}, {0, 0, 0}, 1.0);
        std::vector<Body> bodies{central, particle};
        PostNewtonianGravity pn;
        std::vector<Vec3> acc(2, Vec3{});
        pn.compute(bodies, 0.0, acc);
        const double expected = -(mu_sun / (r0 * r0)) *
            (1.0 - 4.0 * mu_sun / (r0 * c * c));
        std::cout << "  radial a = " << acc[1].x << " expected "
                  << expected << std::endl;
        check("radial 1PN acceleration matches the standard gauge formula",
              std::abs(acc[1].x - expected) < 1.0e-12 * std::abs(expected));
    }

    // --- 3. Perihelion precession vs the analytic 1PN formula -------------
    {
        const double a = 2000.0 * rg;
        const double e = 0.3;
        const double r_p = a * (1.0 - e);
        const double v_p = std::sqrt(mu_sun * (1.0 + e) / (a * (1.0 - e)));
        Body central = make_body({0, 0, 0}, {0, 0, 0}, mass_sun);
        Body particle = make_body({r_p, 0, 0}, {0, v_p, 0}, 1.0);
        std::vector<Body> bodies{central, particle};
        PostNewtonianGravity pn;
        std::vector<Vec3> acc(2, Vec3{});
        pn.compute(bodies, 0.0, acc);
        const double orbital_period = 2.0 * pi *
            std::sqrt(a * a * a / mu_sun);
        const double dt = orbital_period / 4000.0;
        const double duration = 5.5 * orbital_period;
        // Verlet back-step so the initial velocity enters the integrator.
        Vec3 previous_pos = particle.state.pos - particle.state.vel * dt +
                            acc[1] * (0.5 * dt * dt);

        std::vector<double> pericenter_times;
        std::vector<double> pericenter_angles;
        double r_prev2 = particle.state.pos.norm();
        double r_prev = r_prev2;
        for (double t = 0.0; t < duration; t += dt) {
            const Vec3 next_pos = particle.state.pos * 2.0 - previous_pos +
                                  acc[1] * (dt * dt);
            previous_pos = particle.state.pos;
            particle.state.pos = next_pos;
            bodies[1].state.pos = next_pos;
            acc = std::vector<Vec3>(2, Vec3{});
            pn.compute(bodies, t, acc);
            const double r = next_pos.norm();
            if (r_prev2 > r_prev && r > r_prev) {
                pericenter_times.push_back(t);
                pericenter_angles.push_back(
                    std::atan2(next_pos.y, next_pos.x));
            }
            r_prev2 = r_prev;
            r_prev = r;
        }
        // Per-orbit precession from successive pericenter passages.
        double measured_per_orbit = 0.0;
        int orbits = 0;
        if (pericenter_angles.size() >= 2) {
            // Each raw inter-pericenter angle sweep is one full orbit plus
            // the precession advance; average the magnitude directly.
            double swept = 0.0;
            for (std::size_t i = 1; i < pericenter_angles.size(); ++i) {
                swept += std::fabs(pericenter_angles[i] -
                                   pericenter_angles[i - 1]);
            }
            orbits = static_cast<int>(pericenter_angles.size()) - 1;
            measured_per_orbit = swept / orbits;
        }
        const double analytic = 6.0 * pi * mu_sun /
                                (a * (1.0 - e * e) * c * c);
        std::cout << "  precession: measured " << measured_per_orbit
                  << " rad/orbit over " << orbits << " orbits, analytic "
                  << analytic << '\n';
        check("1PN perihelion precession matches the analytic formula",
              orbits >= 2 &&
                  std::abs(measured_per_orbit - analytic) < 0.03 * analytic);
    }

    std::cout << (failures == 0 ? "PASS: 1PN gravity checks"
                                : "FAILURES PRESENT")
              << '\n';
    return failures == 0 ? 0 : 1;
}
