#include "solar/montecarlo.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>

namespace solar {

// Perturb a mission plan's parameters and execute
static MissionReport execute_perturbed(const MissionPlan& base,
                                       std::mt19937& rng,
                                       const MCParams& params) {
    MissionPlan plan = base;

    std::normal_distribution<double> norm(0.0, 1.0);

    // Perturb thrust (affects engine performance)
    if (params.perturb_thrust) {
        double factor = 1.0 + norm(rng) * 0.05;
        plan.thrust_N *= factor;
    }

    // Perturb Isp (directly affects delta-v budget)
    if (params.perturb_isp) {
        double factor = 1.0 + norm(rng) * 0.02;
        plan.isp_s *= factor;
    }

    // Perturb dry mass
    if (params.perturb_mass) {
        double factor = 1.0 + norm(rng) * params.mass_uncertainty;
        plan.dry_mass_kg *= factor;
    }

    // Perturb departure date (±7 days) — changes Lambert trajectory
    if (params.perturb_thrust || params.perturb_isp) {
        double date_shift = norm(rng) * 3.0; // ±3 day 1-sigma
        for (auto& ev : plan.events) {
            if (ev.type == EventType::Launch) {
                ev.jd += date_shift;
            }
        }
    }

    return execute_mission(plan);
}

MCResult run_montecarlo(const MissionPlan& plan, const MCParams& params) {
    MCResult result;
    result.n_total = params.n_samples;

    std::mt19937 rng(params.seed);

    double dv_sum = 0, dv_sq_sum = 0;
    double fuel_sum = 0;
    result.dv_min = 1e30;
    result.dv_max = -1e30;

    for (int i = 0; i < params.n_samples; ++i) {
        auto report = execute_perturbed(plan, rng, params);

        double dv = report.total_dv_kms;
        result.dv_samples.push_back(dv);
        result.success_samples.push_back(report.success);

        dv_sum += dv;
        dv_sq_sum += dv * dv;
        result.dv_min = std::min(result.dv_min, dv);
        result.dv_max = std::max(result.dv_max, dv);
        fuel_sum += report.fuel_remaining_kg;

        if (report.success) {
            result.n_success++;
        } else if (report.failure_reason.find("fuel") != std::string::npos ||
                   report.failure_reason.find("Fuel") != std::string::npos ||
                   report.failure_reason.find("Insufficient") != std::string::npos) {
            result.n_fuel_exhausted++;
        } else {
            result.n_other_failure++;
        }
    }

    result.success_rate = static_cast<double>(result.n_success) / result.n_total;
    result.dv_mean = dv_sum / result.n_total;
    double variance = dv_sq_sum / result.n_total - result.dv_mean * result.dv_mean;
    result.dv_stddev = std::sqrt(std::max(0.0, variance));
    result.fuel_remaining_mean = fuel_sum / result.n_total;

    return result;
}

SensitivityResult run_sensitivity(const MissionPlan& plan, int n_samples) {
    SensitivityResult result;

    // Baseline success rate
    MCParams base_params;
    base_params.n_samples = n_samples;
    base_params.seed = 12345;
    base_params.perturb_thrust = false;
    base_params.perturb_isp = false;
    base_params.perturb_mass = false;
    auto baseline = run_montecarlo(plan, base_params);
    double base_rate = baseline.success_rate;

    // Perturb thrust only
    {
        MCParams p = base_params;
        p.perturb_thrust = true;
        auto r = run_montecarlo(plan, p);
        double delta = std::fabs(r.success_rate - base_rate);
        result.entries.push_back({"thrust (±5%)", delta});
    }

    // Perturb Isp only
    {
        MCParams p = base_params;
        p.perturb_isp = true;
        auto r = run_montecarlo(plan, p);
        double delta = std::fabs(r.success_rate - base_rate);
        result.entries.push_back({"isp (±2%)", delta});
    }

    // Perturb mass only
    {
        MCParams p = base_params;
        p.perturb_mass = true;
        p.mass_uncertainty = 0.03;
        auto r = run_montecarlo(plan, p);
        double delta = std::fabs(r.success_rate - base_rate);
        result.entries.push_back({"dry_mass (±3%)", delta});
    }

    // Perturb fuel mass (via a modified plan)
    {
        MCParams p = base_params;
        // Special: manually perturb fuel in the plan
        std::mt19937 rng(p.seed);
        std::normal_distribution<double> norm(0.0, 1.0);
        int succ = 0;
        for (int i = 0; i < n_samples; ++i) {
            MissionPlan mp = plan;
            mp.fuel_mass_kg *= (1.0 + norm(rng) * 0.02);
            auto rep = execute_mission(mp);
            if (rep.success) succ++;
        }
        double rate = static_cast<double>(succ) / n_samples;
        double delta = std::fabs(rate - base_rate);
        result.entries.push_back({"fuel_mass (±2%)", delta});
    }

    // Sort by sensitivity descending
    std::sort(result.entries.begin(), result.entries.end(),
              [](const SensitivityEntry& a, const SensitivityEntry& b) {
                  return a.sensitivity > b.sensitivity;
              });

    return result;
}

void print_mc_report(const MCResult& result, const MCParams& params) {
    std::cout << std::string(50, '=') << "\n";
    std::cout << "  MONTE CARLO REPORT\n";
    std::cout << std::string(50, '=') << "\n";
    std::cout << std::fixed;
    std::cout << "  Samples: " << result.n_total << " (seed=" << params.seed << ")\n";
    std::cout << "  Perturbations:";
    if (params.perturb_thrust) std::cout << " thrust(±5%)";
    if (params.perturb_isp) std::cout << " isp(±2%)";
    if (params.perturb_mass) std::cout << " mass(±" << std::setprecision(0) << params.mass_uncertainty*100 << "%)";
    std::cout << "\n";
    std::cout << std::string(50, '-') << "\n";
    std::cout << std::setprecision(1);
    std::cout << "  Success rate: " << result.success_rate * 100 << "%"
              << " (" << result.n_success << "/" << result.n_total << ")\n";
    std::cout << "  Fuel exhausted: " << result.n_fuel_exhausted << "\n";
    std::cout << std::setprecision(3);
    std::cout << "  Total dv: " << result.dv_mean << " ± " << result.dv_stddev << " km/s"
              << " [" << result.dv_min << ", " << result.dv_max << "]\n";
    std::cout << "  Fuel remaining: " << std::setprecision(0) << result.fuel_remaining_mean << " kg (mean)\n";
    std::cout << std::string(50, '=') << "\n";
}

} // namespace solar
