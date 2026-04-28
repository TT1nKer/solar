#include "solar/trajectory.h"
#include "solar/ephemeris.h"
#include "solar/kepler.h"
#include "solar/constants.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace solar {

double synodic_period(const Body& a, const Body& b, double mu_central) {
    double mu = (mu_central > 0) ? mu_central : constants::MU_SUN;
    double T1 = orbital_period(a.elements.a, mu);
    double T2 = orbital_period(b.elements.a, mu);
    if (std::fabs(T1 - T2) < 1.0) return 1e30; // same orbit
    return std::fabs(T1 * T2 / (T1 - T2));
}

// ============================================================
// Porkchop
// ============================================================

PorkchopGrid compute_porkchop(const PorkchopParams& params,
                              const std::vector<Body>& bodies) {
    double mu = (params.mu_central > 0) ? params.mu_central : constants::MU_SUN;

    const Body* dep = find_body(bodies, params.departure_body);
    const Body* arr = find_body(bodies, params.arrival_body);

    PorkchopGrid grid;
    if (!dep || !arr) return grid;

    int nd = params.n_depart;
    int na = params.n_arrive;

    // Build date arrays
    double dd = (nd > 1) ? (params.jd_depart_end - params.jd_depart_start) / (nd - 1) : 0;
    double da = (na > 1) ? (params.jd_arrive_end - params.jd_arrive_start) / (na - 1) : 0;

    grid.depart_dates.resize(nd);
    grid.arrive_dates.resize(na);
    for (int i = 0; i < nd; ++i) grid.depart_dates[i] = params.jd_depart_start + i * dd;
    for (int j = 0; j < na; ++j) grid.arrive_dates[j] = params.jd_arrive_start + j * da;

    grid.cells.resize(nd, std::vector<PorkchopCell>(na));

    for (int i = 0; i < nd; ++i) {
        double jd_d = grid.depart_dates[i];
        State s1 = compute_position(*dep, jd_d);

        for (int j = 0; j < na; ++j) {
            double jd_a = grid.arrive_dates[j];
            PorkchopCell& cell = grid.cells[i][j];
            cell.jd_depart = jd_d;
            cell.jd_arrive = jd_a;

            double tof_days = jd_a - jd_d;
            if (tof_days < 10.0) {
                cell.valid = false;
                cell.dv_total = -1;
                continue;
            }

            State s2 = compute_position(*arr, jd_a);
            double tof_sec = tof_days * constants::DAY;

            LambertResult lam = solve_lambert(s1.pos, s2.pos, tof_sec, mu, params.prograde);
            if (!lam.converged) {
                cell.valid = false;
                cell.dv_total = -1;
                continue;
            }

            Vec3 dv1_vec = lam.v1 - s1.vel;
            Vec3 dv2_vec = s2.vel - lam.v2;

            cell.dv1 = dv1_vec.norm();
            cell.dv2 = dv2_vec.norm();
            cell.dv_total = cell.dv1 + cell.dv2;
            cell.tof_days = tof_days;
            cell.valid = true;
        }
    }

    return grid;
}

LaunchWindow find_optimal_window(const PorkchopGrid& grid) {
    LaunchWindow best;
    best.dv_total = std::numeric_limits<double>::max();

    for (const auto& row : grid.cells) {
        for (const auto& cell : row) {
            if (cell.valid && cell.dv_total < best.dv_total) {
                best.jd_depart = cell.jd_depart;
                best.jd_arrive = cell.jd_arrive;
                best.dv_total = cell.dv_total;
                best.dv1 = cell.dv1;
                best.dv2 = cell.dv2;
                best.tof_days = cell.tof_days;
            }
        }
    }

    return best;
}

std::vector<LaunchWindow> find_windows_per_synodic(
    const PorkchopGrid& grid, double synodic_period_days)
{
    if (grid.depart_dates.empty()) return {};

    double jd_start = grid.depart_dates.front();
    double jd_end = grid.depart_dates.back();

    std::vector<LaunchWindow> windows;

    double bin_start = jd_start;
    while (bin_start < jd_end) {
        double bin_end = bin_start + synodic_period_days;

        LaunchWindow best;
        best.dv_total = std::numeric_limits<double>::max();

        for (const auto& row : grid.cells) {
            for (const auto& cell : row) {
                if (!cell.valid) continue;
                if (cell.jd_depart >= bin_start && cell.jd_depart < bin_end
                    && cell.dv_total < best.dv_total) {
                    best.jd_depart = cell.jd_depart;
                    best.jd_arrive = cell.jd_arrive;
                    best.dv_total = cell.dv_total;
                    best.dv1 = cell.dv1;
                    best.dv2 = cell.dv2;
                    best.tof_days = cell.tof_days;
                }
            }
        }

        if (best.dv_total < std::numeric_limits<double>::max()) {
            windows.push_back(best);
        }

        bin_start = bin_end;
    }

    return windows;
}

// ============================================================
// Gravity Assist
// ============================================================

FlybyResult compute_flyby(const FlybyParams& params) {
    FlybyResult result;
    result.valid = false;

    if (params.rp <= params.radius_body) {
        result.altitude_km = params.rp - params.radius_body;
        return result;
    }

    // Hyperbolic excess velocity (body-centered)
    result.v_inf_in = params.v_in_helio - params.v_body;
    result.v_inf_mag = result.v_inf_in.norm();

    if (result.v_inf_mag < 1e-10) return result;

    // Hyperbolic eccentricity and turn angle
    result.e_hyp = 1.0 + params.rp * result.v_inf_mag * result.v_inf_mag / params.mu_body;
    result.turn_angle = 2.0 * std::asin(1.0 / result.e_hyp);
    result.altitude_km = params.rp - params.radius_body;

    // Rotate v_inf_in by turn_angle in the orbital plane
    // The rotation axis is normal to the plane containing v_inf_in and the
    // body's orbital angular momentum (prograde flyby = trailing side)
    Vec3 v_hat = result.v_inf_in.normalized();

    // Use the body velocity direction to define the flyby plane
    Vec3 h_dir = params.v_body.cross(result.v_inf_in);
    double h_mag = h_dir.norm();

    if (h_mag < 1e-10) {
        // v_inf parallel to v_body — degenerate, use arbitrary perpendicular
        Vec3 perp = {0, 0, 1};
        if (std::fabs(v_hat.dot(perp)) > 0.9) perp = {0, 1, 0};
        h_dir = v_hat.cross(perp);
        h_mag = h_dir.norm();
    }

    Vec3 axis = h_dir / h_mag; // rotation axis (normal to flyby plane)

    // Rodrigues rotation: rotate v_inf_in by turn_angle around axis
    double cos_t = std::cos(result.turn_angle);
    double sin_t = std::sin(result.turn_angle);
    Vec3 v_in = result.v_inf_in;
    result.v_inf_out = v_in * cos_t
                     + axis.cross(v_in) * sin_t
                     + axis * axis.dot(v_in) * (1.0 - cos_t);

    // Convert back to heliocentric
    result.v_out_helio = result.v_inf_out + params.v_body;

    // Equivalent free delta-v
    result.dv_equiv = (result.v_out_helio - params.v_in_helio).norm();

    result.valid = true;
    return result;
}

// ============================================================
// Multi-flyby
// ============================================================

MultiFlybyResult evaluate_multi_flyby(const MultiFlybyParams& params,
                                      const std::vector<Body>& bodies) {
    double mu = (params.mu_central > 0) ? params.mu_central : constants::MU_SUN;

    MultiFlybyResult result;
    result.valid = false;
    result.dv_total = 0;

    // Build waypoint list: [departure, flyby1, flyby2, ..., arrival]
    struct Waypoint { std::string name; double jd; };
    std::vector<Waypoint> waypoints;
    waypoints.push_back({params.departure_body, params.jd_depart});
    for (const auto& fb : params.flybys) {
        waypoints.push_back({fb.body_name, fb.jd});
    }
    waypoints.push_back({params.arrival_body, params.jd_arrive});

    if (waypoints.size() < 2) return result;

    // Solve Lambert for each leg
    std::vector<LambertResult> lamberts;
    std::vector<State> body_states;

    for (size_t k = 0; k < waypoints.size(); ++k) {
        const Body* b = find_body(bodies, waypoints[k].name);
        if (!b) return result;
        body_states.push_back(compute_position(*b, waypoints[k].jd));
    }

    for (size_t k = 0; k < waypoints.size() - 1; ++k) {
        double tof = (waypoints[k+1].jd - waypoints[k].jd) * constants::DAY;
        if (tof <= 0) return result;

        LambertResult lam = solve_lambert(
            body_states[k].pos, body_states[k+1].pos, tof, mu, true);
        if (!lam.converged) return result;
        lamberts.push_back(lam);
    }

    // Departure delta-v
    Vec3 dv_dep = lamberts[0].v1 - body_states[0].vel;
    result.dv_departure = dv_dep.norm();

    // Arrival delta-v
    Vec3 dv_arr = body_states.back().vel - lamberts.back().v2;
    result.dv_arrival = dv_arr.norm();

    result.dv_total = result.dv_departure + result.dv_arrival;

    // Process each leg
    for (size_t k = 0; k < lamberts.size(); ++k) {
        MultiFlybyLegResult leg;
        leg.from_body = waypoints[k].name;
        leg.to_body = waypoints[k+1].name;
        leg.jd_depart = waypoints[k].jd;
        leg.jd_arrive = waypoints[k+1].jd;
        leg.tof_days = waypoints[k+1].jd - waypoints[k].jd;
        leg.dv = 0;

        // At flyby nodes (not first or last), compute v_inf mismatch
        if (k > 0) {
            // Incoming v_inf from previous leg
            Vec3 v_in = lamberts[k-1].v2;   // arrival velocity of previous leg
            Vec3 v_out = lamberts[k].v1;     // departure velocity of this leg
            Vec3 v_body = body_states[k].vel;

            Vec3 v_inf_in = v_in - v_body;
            Vec3 v_inf_out = v_out - v_body;

            double v_inf_in_mag = v_inf_in.norm();
            double v_inf_out_mag = v_inf_out.norm();

            // Delta-v deficit: difference in v_inf magnitudes
            // (a perfect unpowered flyby preserves |v_inf|)
            double deficit = std::fabs(v_inf_out_mag - v_inf_in_mag);
            leg.dv = deficit;
            result.dv_total += deficit;
        }

        result.legs.push_back(leg);
    }

    result.total_tof_days = params.jd_arrive - params.jd_depart;
    result.valid = true;
    return result;
}

} // namespace solar
