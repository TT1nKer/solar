#include "solar/nbody.h"
#include "solar/gravity.h"
#include "solar/integrator.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace solar {

void NBodySim::init(std::vector<Body> bodies_in, double t0) {
    bodies = std::move(bodies_in);
    time = t0;
    initialized = false;

    // If no forces were added before init, use default (backward compat)
    if (forces_.empty()) {
        add_default_forces();
    }
}

void NBodySim::add_force(std::unique_ptr<ForceModel> model) {
    if (!model) throw std::invalid_argument("NBodySim::add_force: model must not be null");
    forces_.push_back(std::move(model));
    initialized = false;
}

void NBodySim::remove_force(const std::string& name) {
    forces_.erase(
        std::remove_if(forces_.begin(), forces_.end(),
            [&name](const std::unique_ptr<ForceModel>& f) {
                return f->name() == name;
            }),
        forces_.end());
    initialized = false;
}

void NBodySim::clear_forces() {
    forces_.clear();
    initialized = false;
}

void NBodySim::add_default_forces() {
    add_force(std::make_unique<NewtonianGravity>());
}

std::vector<std::string> NBodySim::force_names() const {
    std::vector<std::string> names;
    for (const auto& f : forces_) {
        names.push_back(f->name());
    }
    return names;
}

std::vector<Vec3> NBodySim::compute_accelerations() const {
    size_t n = bodies.size();
    std::vector<Vec3> acc(n, {0, 0, 0});

    for (const auto& force : forces_) {
        force->compute(bodies, time, acc);
    }

    return acc;
}

void NBodySim::step(double dt) {
    if (!std::isfinite(dt) || dt <= 0.0)
        throw std::invalid_argument("NBodySim::step: dt must be finite and positive");

    if (integrator == IntegratorType::Verlet) {
        for (const auto& force : forces_) {
            if (force->depends_on_velocity()) {
                throw std::logic_error(
                    "NBodySim::step: velocity-dependent force '" + force->name() +
                    "' requires RK4 or DOPRI5");
            }
        }
        if (!initialized) {
            accels = compute_accelerations();
            initialized = true;
        }
    }

    size_t n = bodies.size();

    // Extract data for integrator
    std::vector<State> states(n);
    std::vector<double> masses(n);
    for (size_t i = 0; i < n; ++i) {
        states[i] = bodies[i].state;
        masses[i] = bodies[i].mass;
    }
    double integration_time = time;

    auto accel_func = [this](double stage_time,
                             const std::vector<State>& stage_states,
                             const std::vector<double>& /*masses*/) -> std::vector<Vec3> {
        if (stage_states.size() != this->bodies.size())
            throw std::invalid_argument("NBodySim: integrator stage size mismatch");

        std::vector<Body> stage_bodies = this->bodies;
        for (size_t i = 0; i < stage_states.size(); ++i) {
            stage_bodies[i].state = stage_states[i];
        }

        std::vector<Vec3> acc(stage_states.size(), {0, 0, 0});
        for (const auto& force : this->forces_) {
            force->compute(stage_bodies, stage_time, acc);
        }
        return acc;
    };

    // Select integrator
    switch (integrator) {
        case IntegratorType::DOPRI5: {
            double remaining = dt;
            double h = dt_adaptive > 0.0 ? std::min(dt_adaptive, remaining) : remaining;
            int attempts = 0;
            constexpr int MAX_ATTEMPTS = 100000;

            while (remaining > 0.0) {
                h = std::min(h, remaining);
                auto result = dopri5_step(
                    states, masses, accel_func, integration_time, h, atol, rtol);
                total_steps++;
                attempts++;

                if (!std::isfinite(result.dt_next) || result.dt_next <= 0.0 ||
                    attempts > MAX_ATTEMPTS) {
                    throw std::runtime_error("NBodySim::step: adaptive integrator made no progress");
                }

                if (result.accepted) {
                    states = std::move(result.states);
                    integration_time += h;
                    remaining -= h;
                } else {
                    rejected_steps++;
                }
                h = result.dt_next;
            }
            dt_adaptive = h;
            initialized = false;
            break;
        }
        case IntegratorType::RK4:
            states = rk4_step(states, masses, accel_func, integration_time, dt);
            integration_time += dt;
            total_steps++;
            initialized = false;
            break;
        case IntegratorType::Verlet:
        default:
            verlet_step(states, masses, accel_func, accels, integration_time, dt);
            integration_time += dt;
            total_steps++;
            break;
    }

    // Write back
    for (size_t i = 0; i < n; ++i) {
        bodies[i].state = states[i];
    }
    time = integration_time;

}

void NBodySim::run(double duration, double dt, double output_interval, OutputCallback cb) {
    if (!std::isfinite(duration) || duration < 0.0 ||
        !std::isfinite(dt) || dt <= 0.0 ||
        !std::isfinite(output_interval) || output_interval <= 0.0) {
        throw std::invalid_argument(
            "NBodySim::run: duration must be non-negative; dt and output interval must be positive");
    }

    double next_output = output_interval;
    double elapsed = 0.0;

    if (cb) cb(time, bodies);
    while (elapsed < duration) {
        double h = std::min(dt, duration - elapsed);
        step(h);
        elapsed += h;

        if (elapsed >= next_output || elapsed >= duration) {
            next_output += output_interval;
            if (cb) cb(time, bodies);
        }
    }
}

double NBodySim::total_energy() const {
    double ke = 0.0;
    size_t n = bodies.size();

    for (size_t i = 0; i < n; ++i) {
        ke += 0.5 * (bodies[i].mu / constants::G) * bodies[i].state.vel.norm_sq();
    }

    // Sum potential energy from all force models
    double pe = 0.0;
    for (const auto& force : forces_) {
        pe += force->potential_energy(bodies, time);
    }

    return ke + pe;
}

Vec3 NBodySim::total_angular_momentum() const {
    Vec3 L = {0, 0, 0};
    for (const auto& b : bodies) {
        L += b.state.pos.cross(b.state.vel) * (b.mu / constants::G);
    }
    return L;
}

} // namespace solar
