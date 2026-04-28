#pragma once

#include "solar/body.h"
#include "solar/force_model.h"
#include <vector>
#include <functional>
#include <memory>

namespace solar {

enum class IntegratorType { Verlet, RK4, DOPRI5 };

struct NBodySim {
    std::vector<Body> bodies;
    double time = 0.0;   // current time in seconds from epoch

    // Cached accelerations for Verlet integrator
    std::vector<Vec3> accels;
    bool initialized = false;

    // Integrator selection and adaptive control
    IntegratorType integrator = IntegratorType::Verlet;
    double atol = 1e-10;        // absolute tolerance for adaptive
    double rtol = 1e-10;        // relative tolerance for adaptive
    double dt_adaptive = 0.0;   // current adaptive step size (0 = use initial dt)
    int total_steps = 0;        // step counter for diagnostics
    int rejected_steps = 0;     // rejected step counter

    // Composable force models
    std::vector<std::unique_ptr<ForceModel>> forces_;

    // Initialize simulation from a set of bodies
    // If no forces have been added, automatically adds NewtonianGravity
    void init(std::vector<Body> bodies_in, double t0 = 0.0);

    // Force model management
    void add_force(std::unique_ptr<ForceModel> model);
    void remove_force(const std::string& name);
    void clear_forces();
    void add_default_forces(); // adds NewtonianGravity

    // Compute acceleration on each body by summing all force models
    std::vector<Vec3> compute_accelerations() const;

    // Advance simulation by dt seconds (Verlet integrator)
    void step(double dt);

    // Run for duration seconds, calling callback every output_interval
    using OutputCallback = std::function<void(double t, const std::vector<Body>& bodies)>;
    void run(double duration, double dt, double output_interval, OutputCallback cb);

    // Total mechanical energy (kinetic + potential from all force models)
    double total_energy() const;

    // Total angular momentum vector
    Vec3 total_angular_momentum() const;

    // List active force model names (for logging)
    std::vector<std::string> force_names() const;
};

} // namespace solar
