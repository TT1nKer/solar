#include "solar/gr_correction.h"
#include "solar/gravity.h"
#include "solar/integrator.h"
#include "solar/nbody.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace solar;

static int passed = 0;
static int failed = 0;

static void check(const char* name, bool condition) {
    std::cout << (condition ? "  PASS: " : "  FAIL: ") << name << "\n";
    condition ? passed++ : failed++;
}

int main() {
    std::cout << "=== Dynamics Contract Tests ===\n";

    {
        std::vector<State> states(1);
        std::vector<double> masses = {1.0};
        auto acceleration = [](double time, const std::vector<State>& stage_states,
                               const std::vector<double>&) {
            // a=t verifies stage time; -v verifies stage velocity.
            return std::vector<Vec3>{{time - stage_states[0].vel.x, 0.0, 0.0}};
        };

        double time = 0.0;
        double dt = 0.01;
        for (int i = 0; i < 100; ++i) {
            states = rk4_step(states, masses, acceleration, time, dt);
            time += dt;
        }

        // v' + v = t, v(0)=0; x is its integral.
        double expected_velocity = std::exp(-1.0);
        double expected_position = 0.5 - std::exp(-1.0);
        check("RK4 uses stage time and velocity",
              std::fabs(states[0].vel.x - expected_velocity) < 1e-9 &&
              std::fabs(states[0].pos.x - expected_position) < 1e-9);

        std::vector<State> initial(1);
        auto adaptive = dopri5_step(
            initial, masses, acceleration, 0.0, 0.01, 1e-12, 1e-12);
        double expected_adaptive_velocity = 0.01 - 1.0 + std::exp(-0.01);
        check("DOPRI5 uses stage time and velocity",
              adaptive.accepted &&
              std::fabs(adaptive.states[0].vel.x - expected_adaptive_velocity) < 1e-13);
    }

    {
        Body central;
        central.mu = 1.0;
        Body orbiter;
        orbiter.mu = 1e-12;
        orbiter.state = {{10.0, 0.0, 0.0}, {0.0, 0.3, 0.0}};

        NBodySim simulation;
        simulation.add_force(std::make_unique<NewtonianGravity>());
        simulation.add_force(std::make_unique<GRCorrection>(0));
        simulation.init({central, orbiter});

        bool rejected_incompatible_integrator = false;
        try {
            simulation.step(0.1);
        } catch (const std::logic_error&) {
            rejected_incompatible_integrator = true;
        }
        check("Verlet rejects velocity-dependent forces", rejected_incompatible_integrator);
    }

    {
        NBodySim simulation;
        simulation.integrator = IntegratorType::RK4;
        simulation.init({});
        simulation.run(1.0, 0.3, 1.0, {});
        check("run stops exactly at requested duration", std::fabs(simulation.time - 1.0) < 1e-15);
    }

    std::cout << "=== Results: " << passed << " passed, " << failed << " failed ===\n";
    return failed == 0 ? 0 : 1;
}
