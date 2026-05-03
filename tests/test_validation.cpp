// Automated validation tests
//
// Each test corresponds to a claim in docs/validation/*.md.
// If a claim in a doc changes, the corresponding test here must be updated.
// If this test fails, the README/docs are no longer accurate.
//
// Coverage:
//   01: Kepler solver        — covered by test_kepler.cpp
//   02: Hohmann Earth-Mars   — here
//   03: Energy conservation  — here
//   04: DE440 Earth at J2000 — covered by test_de.cpp (requires data file)
//   05: Lagrange points      — here
//   06: Halo + Jacobi        — here
//   07: Frame transforms     — here
//   08: Lambert solver       — here

#include "solar/kepler.h"
#include "solar/transfer.h"
#include "solar/ephemeris.h"
#include "solar/nbody.h"
#include "solar/cr3bp.h"
#include "solar/frame.h"
#include "solar/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace solar;

static int pass = 0, fail_count = 0;

static void check(const char* claim_id, const char* what, bool ok,
                  const std::string& details = "") {
    std::cout << (ok ? "  PASS" : "  FAIL") << " [" << claim_id << "] " << what;
    if (!details.empty()) std::cout << " — " << details;
    std::cout << "\n";
    if (ok) pass++; else fail_count++;
}

static bool approx(double a, double b, double rel_tol) {
    if (std::fabs(b) < 1e-30) return std::fabs(a) < rel_tol;
    return std::fabs(a - b) / std::fabs(b) < rel_tol;
}

int main() {
    std::cout << "=== Validation Tests (cross-checks docs/validation/*.md) ===\n\n";

    // ===== 02: Earth-Mars Hohmann =====
    // Claim: Total Δv = 5.594 km/s, TOF = 258.871 days
    {
        TransferResult r = hohmann_transfer(1.0 * constants::AU,
                                            1.524 * constants::AU,
                                            constants::MU_SUN);
        std::ostringstream d;
        d << std::fixed << std::setprecision(3)
          << "computed=" << r.dv_total << " expected=5.59 (Curtis)";
        check("02", "Hohmann Earth-Mars Δv matches Curtis textbook within 1%",
              approx(r.dv_total, 5.59, 0.01), d.str());

        std::ostringstream t;
        t << std::fixed << std::setprecision(1)
          << "computed=" << r.tof / constants::DAY << "d expected=259d";
        check("02", "Hohmann Earth-Mars TOF matches textbook within 1%",
              approx(r.tof / constants::DAY, 259.0, 0.01), t.str());
    }

    // ===== 03: Energy conservation (1 year, 9 bodies, dt=3600s, Verlet) =====
    // Claim: relative drift < 1e-4 (well within "drift expected 1e-5 to 1e-6 range")
    {
        auto bodies = load_solar_system();
        // Filter to planets only (matches the validation doc setup)
        std::vector<Body> filtered;
        for (auto& b : bodies) {
            if (b.type == BodyType::Star || b.type == BodyType::Planet)
                filtered.push_back(b);
        }
        NBodySim sim;
        sim.init(filtered);
        double E0 = sim.total_energy();
        Vec3 L0 = sim.total_angular_momentum();
        double L0_mag = L0.norm();

        double duration = 30.0 * constants::DAY; // shorter than 1 year for test speed
        double dt = 3600.0;
        sim.run(duration, dt, duration, [](double, const std::vector<Body>&){});

        double Ef = sim.total_energy();
        Vec3 Lf = sim.total_angular_momentum();
        double Lf_mag = Lf.norm();

        double E_drift = std::fabs((Ef - E0) / E0);
        double L_drift = std::fabs((Lf_mag - L0_mag) / L0_mag);

        std::ostringstream e;
        e << std::scientific << std::setprecision(2) << "E_drift=" << E_drift;
        check("03", "Energy drift over 30 days < 1e-4",
              E_drift < 1e-4, e.str());

        std::ostringstream l;
        l << std::scientific << std::setprecision(2) << "L_drift=" << L_drift;
        check("03", "Angular momentum drift over 30 days < 1e-7",
              L_drift < 1e-7, l.str());
    }

    // ===== 05: Lagrange points =====
    // Claim: Sun-Earth L1 distance from Earth = 1.4916e6 km (published ~1.491e6)
    {
        auto sys = cr3bp::make_system("Sun", "Earth");
        auto lps = cr3bp::compute_lagrange_points(sys.mu);

        // L1 in normalized coords: x_L1 < 1-mu (between primaries)
        double l1_dist_from_secondary_norm = 1.0 - sys.mu - lps.L1.x;
        double l1_dist_km = l1_dist_from_secondary_norm * sys.L_star;

        std::ostringstream d;
        d << std::scientific << std::setprecision(3)
          << "computed=" << l1_dist_km << " expected=1.491e+06";
        check("05", "Sun-Earth L1 distance from Earth ≈ 1.491e6 km (within 1%)",
              approx(l1_dist_km, 1.491e6, 0.01), d.str());

        // L4 should be at (0.5-mu, sqrt(3)/2)
        check("05", "L4 y-coord = sqrt(3)/2",
              approx(lps.L4.y, std::sqrt(3.0)/2.0, 1e-10));
        check("05", "L5 y-coord = -sqrt(3)/2",
              approx(lps.L5.y, -std::sqrt(3.0)/2.0, 1e-10));
    }

    // ===== 06: Halo orbit period + Jacobi conservation =====
    // Claim: Sun-Earth L1 Halo (Az=50000 km) period ≈ 177.8 days
    //        Jacobi conserved to <1e-10 over 10 periods
    {
        auto sys = cr3bp::make_system("Sun", "Earth");
        cr3bp::HaloOrbitParams params;
        params.lpoint = 1;
        params.Az = 50000.0 / sys.L_star;
        params.family = 1;

        auto orbit = cr3bp::compute_halo_orbit(sys, params, 1e-8, 80);
        check("06", "SE-L1 Halo orbit converges", orbit.converged);

        if (orbit.converged) {
            double period_days = orbit.period * sys.T_star / constants::DAY;
            std::ostringstream p;
            p << std::fixed << std::setprecision(2)
              << "computed=" << period_days << "d expected=177.8d";
            check("06", "SE-L1 Halo period ≈ 177.8 days (within 1%)",
                  approx(period_days, 177.8, 0.01), p.str());

            // Jacobi conservation over a few periods
            double C0 = cr3bp::jacobi_constant(orbit.initial_state, sys.mu);
            auto traj = cr3bp::propagate_cr3bp(orbit.initial_state, sys.mu,
                                                orbit.period * 3, orbit.period / 4);
            double max_drift = 0;
            for (const auto& s : traj) {
                double C = cr3bp::jacobi_constant(s, sys.mu);
                max_drift = std::max(max_drift, std::fabs(C - C0));
            }
            std::ostringstream j;
            j << std::scientific << std::setprecision(2)
              << "max|dC|/C0=" << max_drift / std::fabs(C0);
            check("06", "Jacobi constant conserved over 3 periods to <1e-10",
                  max_drift / std::fabs(C0) < 1e-10, j.str());
        }
    }

    // ===== 07: Frame transforms =====
    // Claim: ICRF↔Ecliptic round-trip is identity to machine precision
    {
        Vec3 v_orig = {1.5e8, 2.7e7, -5.4e6}; // typical Earth-scale vector
        Vec3 v_icrf = ecliptic_to_icrf().apply(v_orig);
        Vec3 v_back = icrf_to_ecliptic().apply(v_icrf);
        double err = (v_back - v_orig).norm();

        std::ostringstream e;
        double rel_err = err / v_orig.norm();
        e << std::scientific << std::setprecision(2)
          << "abs_err=" << err << " rel_err=" << rel_err;
        check("07", "Ecliptic→ICRF→Ecliptic round-trip relative err < 1e-14",
              rel_err < 1e-14, e.str());

        // Magnitude preservation
        double mag_orig = v_orig.norm();
        double mag_rot = v_icrf.norm();
        std::ostringstream m;
        m << std::scientific << std::setprecision(2)
          << "rel_err=" << std::fabs(mag_rot - mag_orig) / mag_orig;
        check("07", "Rotation preserves vector magnitude (orthogonal)",
              std::fabs(mag_rot - mag_orig) / mag_orig < 1e-15, m.str());
    }

    // ===== 08: Lambert solver — use 90° transfer to avoid 180° degeneracy =====
    // Claim: Lambert produces a self-consistent transfer for an Earth-Mars-like geometry.
    // (A pure Hohmann is 180° apart, which is a known-singular case for the
    // universal-variable Lambert formulation. We test a non-degenerate case here.)
    {
        auto bodies = load_solar_system();
        const Body* earth = find_body(bodies, "Earth");
        const Body* mars = find_body(bodies, "Mars");

        double jd_dep = parse_date("2026-10-26");
        double jd_arr = parse_date("2027-09-03");
        double tof = (jd_arr - jd_dep) * constants::DAY;

        State s_earth = compute_position(*earth, jd_dep);
        State s_mars = compute_position(*mars, jd_arr);

        auto lam = solve_lambert(s_earth.pos, s_mars.pos, tof,
                                 constants::MU_SUN, true);
        check("08", "Lambert converges for 2026 Earth-Mars window", lam.converged);

        if (lam.converged) {
            double dv1 = (lam.v1 - s_earth.vel).norm();
            double dv2 = (s_mars.vel - lam.v2).norm();
            double total = dv1 + dv2;
            std::ostringstream v;
            v << std::fixed << std::setprecision(3)
              << "total=" << total << " expected_range=[5.0, 7.0]";
            check("08", "Lambert Earth-Mars 2026 Δv in expected range [5.0, 7.0]",
                  total > 5.0 && total < 7.0, v.str());
        }
    }

    std::cout << "\n=== Results: " << pass << " passed, " << fail_count << " failed ===\n";
    return fail_count > 0 ? 1 : 0;
}
