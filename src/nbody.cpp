#include "solar/nbody.h"
#include "solar/gravity.h"
#include "solar/integrator.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>

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
    forces_.push_back(std::move(model));
}

void NBodySim::remove_force(const std::string& name) {
    forces_.erase(
        std::remove_if(forces_.begin(), forces_.end(),
            [&name](const std::unique_ptr<ForceModel>& f) {
                return f->name() == name;
            }),
        forces_.end());
}

void NBodySim::clear_forces() {
    forces_.clear();
}

void NBodySim::add_default_forces() {
    forces_.push_back(std::make_unique<NewtonianGravity>());
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
    if (!initialized) {
        accels = compute_accelerations();
        initialized = true;
    }

    size_t n = bodies.size();

    // Extract data for integrator
    std::vector<State> states(n);
    std::vector<double> masses(n);
    for (size_t i = 0; i < n; ++i) {
        states[i] = bodies[i].state;
        masses[i] = bodies[i].mass;
    }

    // Acceleration function: builds temporary body view with updated positions,
    // then delegates to all force models
    auto accel_func = [this](const std::vector<Vec3>& positions,
                             const std::vector<double>& /*masses*/) -> std::vector<Vec3> {
        size_t n = positions.size();

        // Build temporary bodies with updated positions for force models
        std::vector<Body> temp_bodies = this->bodies;
        for (size_t i = 0; i < n; ++i) {
            temp_bodies[i].state.pos = positions[i];
        }

        std::vector<Vec3> acc(n, {0, 0, 0});
        for (const auto& force : this->forces_) {
            force->compute(temp_bodies, this->time, acc);
        }
        return acc;
    };

    // Select integrator
    switch (integrator) {
        case IntegratorType::DOPRI5: {
            auto result = dopri5_step(states, masses, accel_func, dt, atol, rtol);
            if (result.accepted) {
                states = std::move(result.states);
                dt_adaptive = result.dt_next;
            }
            total_steps++;
            if (!result.accepted) rejected_steps++;
            break;
        }
        case IntegratorType::RK4:
            states = rk4_step(states, masses, accel_func, dt);
            total_steps++;
            break;
        case IntegratorType::Verlet:
        default:
            verlet_step(states, masses, accel_func, accels, dt);
            total_steps++;
            break;
    }

    // Write back
    for (size_t i = 0; i < n; ++i) {
        bodies[i].state = states[i];
    }

    time += dt;
}

void NBodySim::run(double duration, double dt, double output_interval, OutputCallback cb) {
    double next_output = 0.0;
    double elapsed = 0.0;

    if (cb) cb(time, bodies);

    if (integrator == IntegratorType::DOPRI5) {
        // Adaptive stepping
        if (dt_adaptive <= 0.0) dt_adaptive = dt;
        while (elapsed < duration) {
            double h = std::min(dt_adaptive, duration - elapsed);

            // Extract and prepare
            size_t n = bodies.size();
            std::vector<State> states(n);
            std::vector<double> masses(n);
            for (size_t i = 0; i < n; ++i) {
                states[i] = bodies[i].state;
                masses[i] = bodies[i].mass;
            }

            auto accel_func = [this](const std::vector<Vec3>& positions,
                                     const std::vector<double>&) -> std::vector<Vec3> {
                size_t n = positions.size();
                std::vector<Body> temp = this->bodies;
                for (size_t i = 0; i < n; ++i) temp[i].state.pos = positions[i];
                std::vector<Vec3> acc(n, {0,0,0});
                for (const auto& f : this->forces_) f->compute(temp, this->time, acc);
                return acc;
            };

            auto result = dopri5_step(states, masses, accel_func, h, atol, rtol);
            total_steps++;

            if (result.accepted) {
                for (size_t i = 0; i < n; ++i) bodies[i].state = result.states[i];
                time += h;
                elapsed += h;
                dt_adaptive = result.dt_next;

                if (elapsed >= next_output + output_interval) {
                    next_output += output_interval;
                    if (cb) cb(time, bodies);
                }
            } else {
                rejected_steps++;
                dt_adaptive = result.dt_next;
            }
        }
        return;
    }

    // Fixed stepping (Verlet or RK4)
    while (elapsed < duration) {
        step(dt);
        elapsed += dt;

        if (elapsed >= next_output + output_interval) {
            next_output += output_interval;
            if (cb) cb(time, bodies);
        }
    }
}

double NBodySim::total_energy() const {
    double ke = 0.0;
    size_t n = bodies.size();

    for (size_t i = 0; i < n; ++i) {
        ke += 0.5 * bodies[i].mass * bodies[i].state.vel.norm_sq();
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
        L += b.state.pos.cross(b.state.vel) * b.mass;
    }
    return L;
}

} // namespace solar
