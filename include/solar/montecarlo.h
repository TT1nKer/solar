#pragma once

#include "solar/mission.h"
#include <vector>
#include <string>

namespace solar {

struct MCParams {
    int n_samples = 100;
    unsigned seed = 42;
    bool perturb_thrust = true;
    bool perturb_isp = true;
    bool perturb_mass = false;
    double mass_uncertainty = 0.01;  // ±1% dry mass
};

struct MCResult {
    int n_total = 0;
    int n_success = 0;
    int n_fuel_exhausted = 0;
    int n_other_failure = 0;
    double success_rate = 0;
    double dv_mean = 0;
    double dv_stddev = 0;
    double dv_min = 0;
    double dv_max = 0;
    double fuel_remaining_mean = 0;
    std::vector<double> dv_samples;    // total dv per sample
    std::vector<bool> success_samples; // success per sample
};

struct SensitivityEntry {
    std::string param_name;
    double sensitivity;  // |d(success_rate)/d(param_fraction)|
};

struct SensitivityResult {
    std::vector<SensitivityEntry> entries; // sorted by |sensitivity| descending
};

// Run Monte Carlo on a mission plan
MCResult run_montecarlo(const MissionPlan& plan, const MCParams& params);

// Run sensitivity analysis: perturb each parameter individually
SensitivityResult run_sensitivity(const MissionPlan& plan, int n_samples = 50);

// Print MC results
void print_mc_report(const MCResult& result, const MCParams& params);

} // namespace solar
