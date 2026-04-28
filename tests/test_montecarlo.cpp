// TEST-FIRST: Define what "Vehicle + Uncertainty" means before implementing.
//
// This test defines the acceptance criteria:
// 1. A Vehicle has subsystems with nominal + uncertainty parameters
// 2. MonteCarloRunner perturbs parameters, runs mission N times
// 3. Output: success rate, failure breakdown, sensitivity ranking
//
// If this test passes, the capability is real. If not, it's not done.

#include "solar/vehicle.h"
#include "solar/montecarlo.h"
#include "solar/mission.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace solar;

static int pass = 0, fail_count = 0;

static void check(const char* name, bool condition) {
    if (condition) { std::cout << "  PASS: " << name << "\n"; pass++; }
    else { std::cout << "  FAIL: " << name << "\n"; fail_count++; }
}

int main() {
    std::cout << "=== Vehicle + Uncertainty Acceptance Tests ===\n\n";

    // --- Test 1: Vehicle construction with subsystems ---
    std::cout << "Test 1: Vehicle has subsystems with nominal + uncertainty\n";
    {
        Vehicle v;
        v.name = "Mars Cruiser";
        v.dry_mass_kg = 1200;
        v.propulsion.fuel_mass_kg = 8000;
        v.propulsion.thrust_N = 44000;
        v.propulsion.isp_s = 320;
        v.propulsion.thrust_uncertainty = 0.05;  // ±5%
        v.propulsion.isp_uncertainty = 0.02;      // ±2%

        check("vehicle name set", v.name == "Mars Cruiser");
        check("thrust nominal", v.propulsion.thrust_N == 44000);
        check("thrust uncertainty", v.propulsion.thrust_uncertainty == 0.05);
        check("max_dv computable", v.max_delta_v() > 0);

        double dv = v.max_delta_v();
        std::cout << "  max_dv = " << std::fixed << std::setprecision(3) << dv << " km/s\n";
        check("max_dv reasonable (4-8 km/s)", dv > 4.0 && dv < 8.0);
    }

    // --- Test 2: Monte Carlo runner produces statistics ---
    std::cout << "\nTest 2: Monte Carlo runs mission N times with perturbations\n";
    {
        MissionPlan plan = mars_direct_template();
        MCParams params;
        params.n_samples = 100;
        params.seed = 42;
        params.perturb_thrust = true;
        params.perturb_isp = true;

        MCResult result = run_montecarlo(plan, params);

        check("ran 100 samples", result.n_total == 100);
        check("some succeeded", result.n_success > 0);
        check("success_rate in [0,1]", result.success_rate >= 0 && result.success_rate <= 1);
        check("has dv statistics", result.dv_mean > 0);
        check("has dv stddev", result.dv_stddev >= 0);
        check("failure reasons populated", result.n_success + result.n_fuel_exhausted == result.n_total
              || result.n_success + result.n_fuel_exhausted + result.n_other_failure == result.n_total);

        std::cout << std::fixed << std::setprecision(1);
        std::cout << "  Success rate: " << result.success_rate * 100 << "%\n";
        std::cout << "  dv mean: " << std::setprecision(3) << result.dv_mean << " km/s\n";
        std::cout << "  dv stddev: " << result.dv_stddev << " km/s\n";
        std::cout << "  Fuel exhausted: " << result.n_fuel_exhausted << "/" << result.n_total << "\n";
    }

    // --- Test 3: Sensitivity — which parameter matters most ---
    std::cout << "\nTest 3: Sensitivity analysis identifies dominant parameter\n";
    {
        MissionPlan plan = mars_direct_template();
        SensitivityResult sens = run_sensitivity(plan);

        check("has entries", !sens.entries.empty());

        std::cout << "  Parameter sensitivities (|d(success_rate)/d(param)|):\n";
        for (const auto& e : sens.entries) {
            std::cout << "    " << std::left << std::setw(20) << e.param_name
                      << std::right << std::fixed << std::setprecision(4)
                      << " sensitivity=" << e.sensitivity << "\n";
        }

        // The most sensitive parameter should be identifiable
        check("most sensitive identified", !sens.entries.empty() && sens.entries[0].sensitivity > 0);
    }

    std::cout << "\n=== Results: " << pass << " passed, " << fail_count << " failed ===\n";
    return fail_count > 0 ? 1 : 0;
}
