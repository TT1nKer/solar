#include "solar/body.h"
#include "solar/constants.h"
#include "solar/kepler.h"
#include "solar/ephemeris.h"
#include "solar/jpl_de.h"
#include "solar/nbody.h"
#include "solar/transfer.h"
#include "solar/gravity.h"
#include "solar/j2.h"
#include "solar/gr_correction.h"
#include "solar/spherical_harmonics.h"
#include "solar/srp.h"
#include "solar/drag.h"
#include "solar/trajectory.h"
#include "solar/cr3bp.h"
#include "solar/mission.h"
#include "solar/montecarlo.h"
#include "solar/network.h"
#include "solar/vehicle.h"
#include "solar/frame.h"
#include "solar/iau_rotation.h"
#include "solar/time_scale.h"
#include "blackhole_command.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace solar;
using namespace solar::constants;

static const char* body_type_str(BodyType t) {
    switch (t) {
        case BodyType::Star: return "Star";
        case BodyType::Planet: return "Planet";
        case BodyType::Moon: return "Moon";
        case BodyType::Spacecraft: return "Craft";
    }
    return "?";
}

static void print_usage() {
    std::cerr << "Usage: solar <command> [options]\n"
              << "\n"
              << "Commands:\n"
              << "  bodies                                List all bodies with properties\n"
              << "  info <body>                           Detailed info for a body\n"
              << "  moons <planet>                        List moons of a planet\n"
              << "  ephemeris <body> <YYYY-MM-DD>         Position/velocity at date\n"
              << "  orbits                                Orbital elements for planets\n"
              << "  simulate <days> <dt_sec> [interval]   N-body simulation\n"
              << "  transfer <from> <to> <YYYY-MM-DD>     Hohmann transfer\n"
              << "  lambert <from> <to> <depart> <arrive> Lambert problem\n"
              << "  energy <days> <dt_sec>                Energy conservation check\n"
              << "  frame <body> <YYYY-MM-DD> <frame>     Position in given frame\n"
              << "  time <YYYY-MM-DD> [HH:MM]             Time scale conversions\n"
              << "  porkchop <from> <to> <d1> <d2> <a1> <a2> [--grid N]  Porkchop plot (TSV)\n"
              << "  launch-window <from> <to> <d1> <d2> <a1> <a2>        Best launch windows\n"
              << "  flyby <body> <date> <vx> <vy> <vz> <rp_km>          Gravity assist calc\n"
              << "  multi-flyby <from> <to> <dep> <arr> --via <b1> <d1>  Multi-flyby trajectory\n"
              << "  lagrange <primary> <secondary>         Lagrange points L1-L5\n"
              << "  halo <primary> <secondary> <L> <Az_km> Compute Halo orbit\n"
              << "  mission <template>                     Run mission (mars|gateway|grand-tour)\n"
              << "  montecarlo <template> [--samples N]    Monte Carlo uncertainty analysis\n"
              << "  network <YYYY-MM-DD>                   Solar system transfer network\n"
              << "  blackhole <circular|photon> ...        Schwarzschild geodesic simulation\n"
              << "\n"
              << "Physics flags (for simulate/energy):\n"
              << "  --j2             Enable J2 oblateness perturbation\n"
              << "  --harmonics      Enable spherical harmonics (J2-J6)\n"
              << "  --gr             Enable 1PN Schwarzschild correction (requires --adaptive)\n"
              << "  --srp            Enable solar radiation pressure\n"
              << "  --drag           Enable atmospheric drag (requires --adaptive)\n"
              << "  --adaptive       Use DOPRI5 adaptive integrator\n"
              << "  --tol <val>      Adaptive tolerance (default 1e-10)\n"
              << "  --planets-only   Exclude moons from simulation\n";
}

// Helper: check if a flag is present in argv
static bool has_flag(int argc, char* argv[], const char* flag) {
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == flag) return true;
    }
    return false;
}

// Configure force models on NBodySim based on CLI flags
static void configure_forces(NBodySim& sim, int argc, char* argv[],
                             const std::vector<Body>& bodies) {
    bool use_j2 = has_flag(argc, argv, "--j2");
    bool use_gr = has_flag(argc, argv, "--gr");
    bool use_harmonics = has_flag(argc, argv, "--harmonics");
    bool use_srp = has_flag(argc, argv, "--srp");
    bool use_drag = has_flag(argc, argv, "--drag");

    // Always start with gravity
    sim.add_force(std::make_unique<NewtonianGravity>());

    if (use_harmonics) {
        auto sh = SphericalHarmonics::solar_system_defaults(bodies, 6);
        sim.add_force(std::make_unique<SphericalHarmonics>(std::move(sh)));
    } else if (use_j2) {
        auto j2 = J2Perturbation::solar_system_defaults(bodies);
        sim.add_force(std::make_unique<J2Perturbation>(std::move(j2)));
    }

    if (use_gr) {
        sim.add_force(std::make_unique<GRCorrection>(0));
    }

    // SRP and drag apply to spacecraft bodies
    if (use_srp || use_drag) {
        std::vector<SRPTarget> srp_targets;
        std::vector<DragTarget> drag_targets;
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (bodies[i].type == BodyType::Spacecraft) {
                // Default spacecraft properties: 10 m^2 area, Cd=2.2, Cr=1.5
                if (use_srp) srp_targets.push_back({static_cast<int>(i), 1.5, 10.0});
                if (use_drag) drag_targets.push_back({static_cast<int>(i), 2.2, 10.0});
            }
        }
        if (use_srp && !srp_targets.empty())
            sim.add_force(std::make_unique<SolarRadiationPressure>(0, std::move(srp_targets)));
        if (use_drag && !drag_targets.empty())
            sim.add_force(std::make_unique<AtmosphericDrag>(std::move(drag_targets)));
    }
}

// Configure integrator based on CLI flags
static void configure_integrator(NBodySim& sim, int argc, char* argv[]) {
    if (has_flag(argc, argv, "--adaptive")) {
        sim.integrator = IntegratorType::DOPRI5;
        // Check for --tol flag
        for (int i = 2; i < argc - 1; ++i) {
            if (std::string(argv[i]) == "--tol") {
                double tol = std::atof(argv[i + 1]);
                if (tol > 0) { sim.atol = tol; sim.rtol = tol; }
            }
        }
    }
}

// ============================================================
// New Phase 2 commands
// ============================================================

static int cmd_bodies(int /*argc*/, char* /*argv*/[]) {
    auto bodies = load_solar_system();

    std::cout << std::left << std::setw(12) << "Name"
              << std::setw(8) << "Type"
              << std::setw(12) << "Parent"
              << std::right
              << std::setw(12) << "Radius(km)"
              << std::setw(12) << "g(m/s^2)"
              << std::setw(12) << "Vesc(km/s)"
              << std::setw(8) << "Atmo"
              << std::setw(12) << "P(atm)"
              << "\n";
    std::cout << std::string(88, '-') << "\n";

    for (size_t i = 0; i < bodies.size(); ++i) {
        const auto& b = bodies[i];
        std::string parent = (b.parent_index >= 0) ? bodies[b.parent_index].name : "-";

        std::cout << std::left << std::setw(12) << b.name
                  << std::setw(8) << body_type_str(b.type)
                  << std::setw(12) << parent
                  << std::right << std::fixed
                  << std::setw(12) << std::setprecision(1) << b.properties.radius_km
                  << std::setw(12) << std::setprecision(2) << b.properties.surface_gravity
                  << std::setw(12) << std::setprecision(3) << b.properties.escape_velocity
                  << std::setw(8) << (b.properties.has_atmosphere ? "yes" : "no")
                  << std::setw(12) << std::setprecision(3) << b.properties.atm_surface_pressure
                  << "\n";
    }

    return 0;
}

static int cmd_info(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: solar info <body>\n";
        return 1;
    }

    auto bodies = load_solar_system();
    const Body* body = find_body(bodies, argv[2]);
    if (!body) {
        std::cerr << "Unknown body: " << argv[2] << "\n";
        return 1;
    }

    std::string parent = (body->parent_index >= 0) ? bodies[body->parent_index].name : "none";

    std::cout << "=== " << body->name << " ===\n";
    std::cout << "Type:              " << body_type_str(body->type) << "\n";
    std::cout << "Parent:            " << parent << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Mass:              " << body->mass << " kg\n";
    std::cout << "GM (mu):           " << body->mu << " km^3/s^2\n";

    std::cout << "\nPhysical Properties:\n";
    std::cout << std::fixed;
    std::cout << "  Radius:          " << std::setprecision(1) << body->properties.radius_km << " km\n";
    std::cout << "  Surface gravity: " << std::setprecision(2) << body->properties.surface_gravity << " m/s^2\n";
    std::cout << "  Escape velocity: " << std::setprecision(3) << body->properties.escape_velocity << " km/s\n";
    std::cout << "  Atmosphere:      " << (body->properties.has_atmosphere ? "yes" : "no") << "\n";
    if (body->properties.has_atmosphere) {
        std::cout << "  Surface pressure:" << std::setprecision(3) << body->properties.atm_surface_pressure << " atm\n";
    }

    if (body->type != BodyType::Star) {
        double mu_p = (body->mu_parent > 0.0) ? body->mu_parent : MU_SUN;
        std::cout << "\nOrbital Elements (relative to " << parent << "):\n";
        std::cout << std::fixed;
        if (body->type == BodyType::Planet) {
            std::cout << "  Semi-major axis: " << std::setprecision(4) << body->elements.a / AU << " AU ("
                      << std::scientific << body->elements.a << " km)\n";
        } else {
            std::cout << "  Semi-major axis: " << std::scientific << std::setprecision(4) << body->elements.a << " km\n";
        }
        std::cout << std::fixed;
        std::cout << "  Eccentricity:    " << std::setprecision(6) << body->elements.e << "\n";
        std::cout << "  Inclination:     " << std::setprecision(4) << body->elements.i * RAD2DEG << " deg\n";
        std::cout << "  Orbital period:  " << std::setprecision(2) << orbital_period(body->elements.a, mu_p) / DAY << " days\n";
    }

    // List children (moons)
    std::vector<std::string> children;
    int body_idx = -1;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (&bodies[i] == body) body_idx = static_cast<int>(i);
    }
    if (body_idx >= 0) {
        for (const auto& b : bodies) {
            if (b.parent_index == body_idx) children.push_back(b.name);
        }
    }
    if (!children.empty()) {
        std::cout << "\nMoons/Children:\n";
        for (const auto& c : children) {
            std::cout << "  - " << c << "\n";
        }
    }

    return 0;
}

static int cmd_moons(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: solar moons <planet>\n";
        return 1;
    }

    auto bodies = load_solar_system();
    const Body* planet = find_body(bodies, argv[2]);
    if (!planet) {
        std::cerr << "Unknown body: " << argv[2] << "\n";
        return 1;
    }

    // Find planet's index
    int planet_idx = -1;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (&bodies[i] == planet) { planet_idx = static_cast<int>(i); break; }
    }

    std::cout << "Moons of " << planet->name << ":\n\n";
    std::cout << std::left << std::setw(12) << "Name"
              << std::right
              << std::setw(14) << "a (km)"
              << std::setw(12) << "e"
              << std::setw(14) << "Period (d)"
              << std::setw(12) << "R (km)"
              << std::setw(12) << "g (m/s^2)"
              << "\n";
    std::cout << std::string(76, '-') << "\n";

    bool found = false;
    for (const auto& b : bodies) {
        if (b.parent_index != planet_idx || b.type != BodyType::Moon) continue;
        found = true;

        double period_d = orbital_period(b.elements.a, planet->mu) / DAY;

        std::cout << std::left << std::setw(12) << b.name
                  << std::right << std::fixed
                  << std::setw(14) << std::setprecision(0) << b.elements.a
                  << std::setw(12) << std::setprecision(4) << b.elements.e
                  << std::setw(14) << std::setprecision(2) << period_d
                  << std::setw(12) << std::setprecision(1) << b.properties.radius_km
                  << std::setw(12) << std::setprecision(4) << b.properties.surface_gravity
                  << "\n";
    }

    if (!found) {
        std::cout << "(no moons loaded for " << planet->name << ")\n";
    }

    return 0;
}

// ============================================================
// Updated Phase 1 commands
// ============================================================

static int cmd_ephemeris(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: solar ephemeris <body> <YYYY-MM-DD>\n";
        return 1;
    }

    std::string body_name = argv[2];
    double jd = parse_date(argv[3]);

    auto bodies = load_solar_system();
    const Body* body = find_body(bodies, body_name);
    if (!body) {
        std::cerr << "Unknown body: " << body_name << "\n";
        return 1;
    }

    // For moons, show both parent-relative and heliocentric
    bool is_moon = (body->type == BodyType::Moon);
    State s = compute_position(*body, jd);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Body: " << body->name << "\n";
    std::cout << "Date: " << argv[3] << " (JD " << std::fixed << std::setprecision(1) << jd << ")\n";

    if (is_moon) {
        std::string parent_name = bodies[body->parent_index].name;
        std::cout << "\nRelative to " << parent_name << ":\n";
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "  Position (km):   " << std::setw(14) << s.pos.x
                  << "  " << std::setw(14) << s.pos.y
                  << "  " << std::setw(14) << s.pos.z << "\n";
        std::cout << "  Velocity (km/s): " << std::setw(14) << s.vel.x
                  << "  " << std::setw(14) << s.vel.y
                  << "  " << std::setw(14) << s.vel.z << "\n";
        std::cout << "  Distance: " << s.pos.norm() << " km\n";
        std::cout << "  Speed:    " << std::fixed << std::setprecision(3) << s.vel.norm() << " km/s\n";

        State sh = compute_heliocentric(*body, bodies, jd);
        std::cout << "\nHeliocentric:\n";
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "  Position (km):   " << std::setw(14) << sh.pos.x
                  << "  " << std::setw(14) << sh.pos.y
                  << "  " << std::setw(14) << sh.pos.z << "\n";
        std::cout << "  Velocity (km/s): " << std::setw(14) << sh.vel.x
                  << "  " << std::setw(14) << sh.vel.y
                  << "  " << std::setw(14) << sh.vel.z << "\n";
        double dist = sh.pos.norm();
        std::cout << "  Distance from Sun: " << dist << " km ("
                  << std::fixed << std::setprecision(4) << dist / AU << " AU)\n";
    } else {
        double dist = s.pos.norm();
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "Position (km):   " << std::setw(14) << s.pos.x
                  << "  " << std::setw(14) << s.pos.y
                  << "  " << std::setw(14) << s.pos.z << "\n";
        std::cout << "Velocity (km/s): " << std::setw(14) << s.vel.x
                  << "  " << std::setw(14) << s.vel.y
                  << "  " << std::setw(14) << s.vel.z << "\n";
        std::cout << "Distance from Sun: " << std::setprecision(4) << dist
                  << " km (" << std::fixed << std::setprecision(4) << dist / AU << " AU)\n";
        std::cout << "Speed: " << std::fixed << std::setprecision(3) << s.vel.norm() << " km/s\n";
    }

    return 0;
}

static int cmd_orbits(int /*argc*/, char* /*argv*/[]) {
    auto bodies = load_solar_system();

    std::cout << std::left << std::setw(10) << "Body"
              << std::right
              << std::setw(14) << "a (AU)"
              << std::setw(12) << "e"
              << std::setw(12) << "i (deg)"
              << std::setw(14) << "Omega (deg)"
              << std::setw(14) << "omega (deg)"
              << std::setw(14) << "Period (yr)"
              << "\n";
    std::cout << std::string(90, '-') << "\n";

    for (const auto& b : bodies) {
        if (b.type != BodyType::Planet) continue;

        double a_au = b.elements.a / AU;
        double period_yr = orbital_period(b.elements.a, MU_SUN) / YEAR;

        std::cout << std::left << std::setw(10) << b.name
                  << std::right << std::fixed
                  << std::setw(14) << std::setprecision(6) << a_au
                  << std::setw(12) << std::setprecision(6) << b.elements.e
                  << std::setw(12) << std::setprecision(4) << b.elements.i * RAD2DEG
                  << std::setw(14) << std::setprecision(4) << b.elements.Omega * RAD2DEG
                  << std::setw(14) << std::setprecision(4) << b.elements.omega * RAD2DEG
                  << std::setw(14) << std::setprecision(3) << period_yr
                  << "\n";
    }

    return 0;
}

static int cmd_simulate(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: solar simulate <days> <dt_sec> [output_interval_days] [--planets-only]\n";
        return 1;
    }

    double days = std::atof(argv[2]);
    double dt = std::atof(argv[3]);
    double output_days = (argc >= 5 && argv[4][0] != '-') ? std::atof(argv[4]) : days / 10.0;

    bool planets_only = false;
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--planets-only") planets_only = true;
    }

    double duration = days * DAY;
    double output_interval = output_days * DAY;

    auto bodies = load_solar_system();

    // For N-body sim, convert moon states to heliocentric
    for (auto& b : bodies) {
        if (b.type == BodyType::Moon) {
            // Convert parent-relative state to heliocentric
            State parent_s = compute_position(bodies[b.parent_index], J2000);
            b.state.pos = b.state.pos + parent_s.pos;
            b.state.vel = b.state.vel + parent_s.vel;
        }
    }

    // Filter to planets only if requested
    if (planets_only) {
        std::vector<Body> filtered;
        for (auto& b : bodies) {
            if (b.type == BodyType::Star || b.type == BodyType::Planet) {
                filtered.push_back(b);
            }
        }
        bodies = filtered;
    }

    NBodySim sim;
    configure_forces(sim, argc, argv, bodies);
    configure_integrator(sim, argc, argv);
    sim.init(bodies);

    // List active forces
    auto fnames = sim.force_names();
    std::cout << "# N-body simulation: " << days << " days, dt=" << dt << "s"
              << ", output every " << output_days << " days"
              << ", " << bodies.size() << " bodies\n";
    std::cout << "# Forces:";
    for (const auto& fn : fnames) std::cout << " " << fn;
    std::cout << "\n";
    std::cout << "# time(days)  body  x(km)  y(km)  z(km)  vx(km/s)  vy(km/s)  vz(km/s)\n";

    sim.run(duration, dt, output_interval,
        [](double t, const std::vector<Body>& bodies) {
            double t_days = t / DAY;
            for (const auto& b : bodies) {
                std::cout << std::fixed << std::setprecision(2) << t_days
                          << "\t" << b.name
                          << std::scientific << std::setprecision(6)
                          << "\t" << b.state.pos.x
                          << "\t" << b.state.pos.y
                          << "\t" << b.state.pos.z
                          << "\t" << b.state.vel.x
                          << "\t" << b.state.vel.y
                          << "\t" << b.state.vel.z
                          << "\n";
            }
        });

    return 0;
}

static int cmd_transfer(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: solar transfer <from> <to> <YYYY-MM-DD>\n";
        return 1;
    }

    std::string from_name = argv[2];
    std::string to_name = argv[3];

    auto bodies = load_solar_system();
    const Body* from = find_body(bodies, from_name);
    const Body* to = find_body(bodies, to_name);

    if (!from) { std::cerr << "Unknown body: " << from_name << "\n"; return 1; }
    if (!to) { std::cerr << "Unknown body: " << to_name << "\n"; return 1; }

    // Use mean orbital radii for Hohmann calculation
    double r1 = from->elements.a;
    double r2 = to->elements.a;

    // For moons, use parent planet's orbit
    if (from->type == BodyType::Moon) r1 = bodies[from->parent_index].elements.a;
    if (to->type == BodyType::Moon) r2 = bodies[to->parent_index].elements.a;

    TransferResult result = hohmann_transfer(r1, r2, MU_SUN);

    std::cout << "Hohmann Transfer: " << from->name << " -> " << to->name << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Departure orbit: " << r1 << " km (" << std::fixed << std::setprecision(4) << r1 / AU << " AU)\n";
    std::cout << "Arrival orbit:   " << std::scientific << std::setprecision(4) << r2 << " km (" << std::fixed << std::setprecision(4) << r2 / AU << " AU)\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Delta-v1 (departure): " << result.dv1 << " km/s\n";
    std::cout << "Delta-v2 (arrival):   " << result.dv2 << " km/s\n";
    std::cout << "Total delta-v:        " << result.dv_total << " km/s\n";
    std::cout << "Time of flight:       " << result.tof / DAY << " days\n";
    std::cout << "Transfer semi-major:  " << std::scientific << std::setprecision(4) << result.a_transfer << " km\n";
    std::cout << "Transfer eccentricity: " << std::fixed << std::setprecision(6) << result.e_transfer << "\n";

    return 0;
}

static int cmd_lambert(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: solar lambert <from> <to> <depart-YYYY-MM-DD> <arrive-YYYY-MM-DD>\n";
        return 1;
    }

    std::string from_name = argv[2];
    std::string to_name = argv[3];
    double jd_depart = parse_date(argv[4]);
    double jd_arrive = parse_date(argv[5]);

    auto bodies = load_solar_system();
    const Body* from = find_body(bodies, from_name);
    const Body* to = find_body(bodies, to_name);

    if (!from) { std::cerr << "Unknown body: " << from_name << "\n"; return 1; }
    if (!to) { std::cerr << "Unknown body: " << to_name << "\n"; return 1; }

    TransferResult result = plan_transfer(*from, *to, jd_depart, jd_arrive, MU_SUN);

    if (result.dv1 < 0) {
        std::cerr << "Lambert solver did not converge.\n";
        return 1;
    }

    double tof_days = result.tof / DAY;

    std::cout << "Lambert Transfer: " << from->name << " -> " << to->name << "\n";
    std::cout << "Departure: " << argv[4] << " (JD " << std::fixed << std::setprecision(1) << jd_depart << ")\n";
    std::cout << "Arrival:   " << argv[5] << " (JD " << std::fixed << std::setprecision(1) << jd_arrive << ")\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Time of flight:       " << tof_days << " days\n";
    std::cout << "Delta-v1 (departure): " << result.dv1 << " km/s\n";
    std::cout << "Delta-v2 (arrival):   " << result.dv2 << " km/s\n";
    std::cout << "Total delta-v:        " << result.dv_total << " km/s\n";
    std::cout << "Transfer semi-major:  " << std::scientific << std::setprecision(4) << result.a_transfer << " km\n";
    std::cout << "Transfer eccentricity: " << std::fixed << std::setprecision(6) << result.e_transfer << "\n";

    return 0;
}

static int cmd_energy(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: solar energy <days> <dt_sec> [--j2] [--gr]\n";
        return 1;
    }

    double days = std::atof(argv[2]);
    double dt = std::atof(argv[3]);
    double duration = days * DAY;
    double output_interval = duration / 20.0;

    auto bodies = load_solar_system();

    // Convert moon states to heliocentric for N-body
    for (auto& b : bodies) {
        if (b.type == BodyType::Moon) {
            State parent_s = compute_position(bodies[b.parent_index], J2000);
            b.state.pos = b.state.pos + parent_s.pos;
            b.state.vel = b.state.vel + parent_s.vel;
        }
    }

    if (has_flag(argc, argv, "--planets-only")) {
        std::vector<Body> planets;
        for (const auto& body : bodies) {
            if (body.type == BodyType::Star || body.type == BodyType::Planet)
                planets.push_back(body);
        }
        bodies = std::move(planets);
    }

    NBodySim sim;
    configure_forces(sim, argc, argv, bodies);
    configure_integrator(sim, argc, argv);
    sim.init(bodies);

    auto fnames = sim.force_names();

    double E0 = sim.total_energy();
    Vec3 L0 = sim.total_angular_momentum();
    double L0_mag = L0.norm();

    std::cout << "# Energy conservation test: " << days << " days, dt=" << dt << "s"
              << ", " << bodies.size() << " bodies\n";
    std::cout << "# Forces:";
    for (const auto& fn : fnames) std::cout << " " << fn;
    std::cout << "\n";
    std::cout << "# Initial energy: " << std::scientific << std::setprecision(6) << E0 << "\n";
    std::cout << "# Initial |L|:    " << L0_mag << "\n";

    sim.run(duration, dt, output_interval,
        [](double /*t*/, const std::vector<Body>& /*bodies*/) {});

    double Ef = sim.total_energy();
    Vec3 Lf = sim.total_angular_momentum();
    double Lf_mag = Lf.norm();

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "# Final energy:   " << Ef << "\n";
    std::cout << "# Final |L|:      " << Lf_mag << "\n";
    std::cout << "# Relative energy drift:  " << std::fabs((Ef - E0) / E0) << "\n";
    std::cout << "# Relative |L| drift:     " << std::fabs((Lf_mag - L0_mag) / L0_mag) << "\n";

    std::cout << "# Steps: " << sim.total_steps
              << " (rejected: " << sim.rejected_steps << ")\n";

    double drift = std::fabs((Ef - E0) / E0);
    if (drift < 1e-10) {
        std::cout << "# PASS: energy conserved within 1e-10\n";
    } else if (drift < 1e-6) {
        std::cout << "# PASS: energy conserved within 1e-6\n";
    } else if (drift < 1e-3) {
        std::cout << "# WARN: energy drift " << drift << " (consider smaller dt or --adaptive)\n";
    } else {
        std::cout << "# FAIL: energy drift " << drift << " (reduce dt or check integrator)\n";
    }

    return 0;
}

static int cmd_frame(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: solar frame <body> <YYYY-MM-DD> <frame>\n"
                  << "Frames: ecliptic, icrf, bodyfixed\n";
        return 1;
    }

    std::string body_name = argv[2];
    double jd = parse_date(argv[3]);
    Frame target_frame = parse_frame(argv[4]);

    auto bodies = load_solar_system();
    const Body* body = find_body(bodies, body_name);
    if (!body) { std::cerr << "Unknown body: " << body_name << "\n"; return 1; }

    // Get heliocentric ecliptic state
    State s_ecl;
    if (body->type == BodyType::Moon) {
        s_ecl = compute_heliocentric(*body, bodies, jd);
    } else {
        s_ecl = compute_position(*body, jd);
    }

    std::cout << "Body: " << body->name << "\n";
    std::cout << "Date: " << argv[3] << " (JD " << std::fixed << std::setprecision(1) << jd << ")\n";
    std::cout << "Frame: " << frame_name(target_frame) << "\n\n";

    if (target_frame == Frame::BodyFixed) {
        // Body-fixed: need to express position relative to a central body
        // For planets, show position relative to Sun (not very useful)
        // More useful: show a satellite's position relative to the planet in body-fixed
        // For now, show the body's own rotation state
        std::cout << "Note: BodyFixed shows coordinates relative to body center.\n";
        std::cout << "For a planet, this shows the Sun's position in the planet's frame.\n\n";

        // Reverse: show Sun relative to body
        int de_id = body_name_to_de_id(body->name);
        State s_body_centered = {{-s_ecl.pos.x, -s_ecl.pos.y, -s_ecl.pos.z},
                                 {-s_ecl.vel.x, -s_ecl.vel.y, -s_ecl.vel.z}};
        State s_bf = transform_state(s_body_centered, Frame::EclipticJ2000,
                                     Frame::BodyFixed, jd, de_id);
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "Sun in " << body->name << " body-fixed:\n";
        std::cout << "  Position (km): " << s_bf.pos.x << "  " << s_bf.pos.y << "  " << s_bf.pos.z << "\n";

        LatLonAlt lla = cartesian_to_geodetic(s_bf.pos, body->properties.radius_km);
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  Sub-solar point: lat=" << lla.lat_deg << "° lon=" << lla.lon_deg << "°\n";
    } else {
        State s_out = transform_state(s_ecl, Frame::EclipticJ2000, target_frame, jd);
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "Position (km):   " << std::setw(14) << s_out.pos.x
                  << "  " << std::setw(14) << s_out.pos.y
                  << "  " << std::setw(14) << s_out.pos.z << "\n";
        std::cout << "Velocity (km/s): " << std::setw(14) << s_out.vel.x
                  << "  " << std::setw(14) << s_out.vel.y
                  << "  " << std::setw(14) << s_out.vel.z << "\n";
        std::cout << "Distance: " << s_out.pos.norm() << " km\n";
        std::cout << "Speed:    " << std::fixed << std::setprecision(3) << s_out.vel.norm() << " km/s\n";
    }

    return 0;
}

static int cmd_time(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: solar time <YYYY-MM-DD> [HH:MM]\n";
        return 1;
    }

    double jd_utc = parse_date(argv[2]);
    // Parse optional time
    if (argc >= 4) {
        std::string hm = argv[3];
        int h = 0, m = 0;
        if (sscanf(hm.c_str(), "%d:%d", &h, &m) >= 1) {
            jd_utc += (h + m / 60.0) / 24.0;
        }
    }

    Epoch utc = Epoch::from_jd(jd_utc, TimeScale::UTC);
    Epoch tai = utc.to(TimeScale::TAI);
    Epoch tt = utc.to(TimeScale::TT);
    Epoch tdb = utc.to(TimeScale::TDB);

    int year, month, day;
    double hour;
    jd_to_calendar(jd_utc, year, month, day, hour);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Input: " << year << "-" << std::setfill('0') << std::setw(2) << month
              << "-" << std::setw(2) << day << std::setfill(' ')
              << " " << std::setprecision(2) << hour << "h UTC\n\n";
    std::cout << std::setprecision(6);
    std::cout << "UTC:  JD " << utc.jd << "\n";
    std::cout << "TAI:  JD " << tai.jd << "  (UTC + " << leap_seconds_at(jd_utc) << " leap seconds)\n";
    std::cout << "TT:   JD " << tt.jd << "  (TAI + 32.184s)\n";
    std::cout << "TDB:  JD " << tdb.jd << "  (TT + periodic ~1.7ms)\n";

    double dt_tt_utc = (tt.jd - utc.jd) * constants::DAY;
    std::cout << "\nTT - UTC = " << std::setprecision(3) << dt_tt_utc << " seconds\n";
    double dt_tdb_tt = (tdb.jd - tt.jd) * constants::DAY * 1000.0;
    std::cout << "TDB - TT = " << std::setprecision(4) << dt_tdb_tt << " milliseconds\n";

    return 0;
}

static int cmd_porkchop(int argc, char* argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: solar porkchop <from> <to> <depart-start> <depart-end>"
                  << " <arrive-start> <arrive-end> [--grid N]\n";
        return 1;
    }
    auto bodies = load_solar_system();
    PorkchopParams p;
    p.departure_body = argv[2];
    p.arrival_body = argv[3];
    p.jd_depart_start = parse_date(argv[4]);
    p.jd_depart_end = parse_date(argv[5]);
    p.jd_arrive_start = parse_date(argv[6]);
    p.jd_arrive_end = parse_date(argv[7]);
    p.n_depart = 100; p.n_arrive = 100;
    for (int i = 8; i < argc - 1; ++i) {
        if (std::string(argv[i]) == "--grid") { p.n_depart = p.n_arrive = std::atoi(argv[i+1]); }
    }

    std::cerr << "Computing porkchop: " << p.departure_body << " -> " << p.arrival_body
              << " (" << p.n_depart << "x" << p.n_arrive << " grid)...\n";

    auto grid = compute_porkchop(p, bodies);
    auto opt = find_optimal_window(grid);

    // TSV to stdout
    std::cout << "# porkchop: " << p.departure_body << " -> " << p.arrival_body << "\n";
    std::cout << "# depart_jd\tarrive_jd\tdv_total\tdv1\tdv2\ttof_days\n";
    for (const auto& row : grid.cells) {
        for (const auto& c : row) {
            if (!c.valid) continue;
            std::cout << std::fixed << std::setprecision(1) << c.jd_depart
                      << "\t" << c.jd_arrive
                      << std::setprecision(3) << "\t" << c.dv_total
                      << "\t" << c.dv1 << "\t" << c.dv2
                      << std::setprecision(1) << "\t" << c.tof_days << "\n";
        }
    }

    // Summary to stderr
    int y, m, d; double h;
    jd_to_calendar(opt.jd_depart, y, m, d, h);
    std::cerr << std::fixed << std::setprecision(3);
    std::cerr << "Optimal: depart " << y << "-" << std::setfill('0')
              << std::setw(2) << m << "-" << std::setw(2) << d << std::setfill(' ');
    jd_to_calendar(opt.jd_arrive, y, m, d, h);
    std::cerr << ", arrive " << y << "-" << std::setfill('0')
              << std::setw(2) << m << "-" << std::setw(2) << d << std::setfill(' ')
              << ", dv=" << opt.dv_total << " km/s"
              << ", TOF=" << std::setprecision(0) << opt.tof_days << "d\n";
    return 0;
}

static int cmd_launch_window(int argc, char* argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: solar launch-window <from> <to> <depart-start> <depart-end>"
                  << " <arrive-start> <arrive-end> [--grid N]\n";
        return 1;
    }
    auto bodies = load_solar_system();
    PorkchopParams p;
    p.departure_body = argv[2]; p.arrival_body = argv[3];
    p.jd_depart_start = parse_date(argv[4]); p.jd_depart_end = parse_date(argv[5]);
    p.jd_arrive_start = parse_date(argv[6]); p.jd_arrive_end = parse_date(argv[7]);
    p.n_depart = 200; p.n_arrive = 200;
    for (int i = 8; i < argc - 1; ++i) {
        if (std::string(argv[i]) == "--grid") { p.n_depart = p.n_arrive = std::atoi(argv[i+1]); }
    }

    auto grid = compute_porkchop(p, bodies);
    auto opt = find_optimal_window(grid);

    // Synodic period
    const Body* dep = find_body(bodies, p.departure_body);
    const Body* arr = find_body(bodies, p.arrival_body);
    double syn = synodic_period(*dep, *arr, constants::MU_SUN);
    double syn_days = syn / constants::DAY;
    auto windows = find_windows_per_synodic(grid, syn_days);

    std::cout << "Launch Windows: " << p.departure_body << " -> " << p.arrival_body << "\n";
    std::cout << "Synodic period: " << std::fixed << std::setprecision(0) << syn_days << " days\n\n";
    std::cout << std::left << std::setw(14) << "Depart" << std::setw(14) << "Arrive"
              << std::right << std::setw(10) << "dv(km/s)" << std::setw(10) << "TOF(d)" << "\n";
    std::cout << std::string(48, '-') << "\n";

    for (const auto& w : windows) {
        int y, m, d; double h;
        jd_to_calendar(w.jd_depart, y, m, d, h);
        std::cout << std::right << y << "-" << std::setfill('0') << std::setw(2) << m
                  << "-" << std::setw(2) << d << std::setfill(' ') << "    ";
        jd_to_calendar(w.jd_arrive, y, m, d, h);
        std::cout << y << "-" << std::setfill('0') << std::setw(2) << m
                  << "-" << std::setw(2) << d << std::setfill(' ');
        std::cout << std::right << std::fixed
                  << std::setw(10) << std::setprecision(3) << w.dv_total
                  << std::setw(10) << std::setprecision(0) << w.tof_days << "\n";
    }

    std::cout << "\nGlobal optimum: dv=" << std::setprecision(3) << opt.dv_total << " km/s\n";
    return 0;
}

static int cmd_flyby(int argc, char* argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: solar flyby <body> <date> <vx> <vy> <vz> <rp_km>\n";
        return 1;
    }
    auto bodies = load_solar_system();
    const Body* fb_body = find_body(bodies, argv[2]);
    if (!fb_body) { std::cerr << "Unknown body: " << argv[2] << "\n"; return 1; }

    double jd = parse_date(argv[3]);
    FlybyParams fp;
    fp.mu_body = fb_body->mu;
    fp.radius_body = fb_body->properties.radius_km;
    State fb_state = compute_position(*fb_body, jd);
    fp.v_body = fb_state.vel;
    fp.v_in_helio = {std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
    fp.rp = std::atof(argv[7]);

    auto result = compute_flyby(fp);

    std::cout << "Gravity Assist: " << fb_body->name << "\n";
    std::cout << "Date: " << argv[3] << " (JD " << std::fixed << std::setprecision(1) << jd << ")\n\n";

    if (!result.valid) {
        std::cerr << "Invalid flyby: rp < body radius (alt=" << std::setprecision(0)
                  << result.altitude_km << " km)\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "V_inf (in):     " << result.v_inf_mag << " km/s\n";
    std::cout << "Turn angle:     " << result.turn_angle * constants::RAD2DEG << " deg\n";
    std::cout << "Eccentricity:   " << std::setprecision(4) << result.e_hyp << "\n";
    std::cout << "Periapsis alt:  " << std::setprecision(0) << result.altitude_km << " km\n";
    std::cout << "Equiv delta-v:  " << std::setprecision(3) << result.dv_equiv << " km/s\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "V_out (helio):  " << result.v_out_helio.x << "  "
              << result.v_out_helio.y << "  " << result.v_out_helio.z << " km/s\n";
    return 0;
}

static int cmd_multi_flyby(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: solar multi-flyby <from> <to> <depart> <arrive>"
                  << " --via <body1> <date1> [<body2> <date2> ...]\n";
        return 1;
    }
    auto bodies = load_solar_system();
    MultiFlybyParams p;
    p.departure_body = argv[2]; p.arrival_body = argv[3];
    p.jd_depart = parse_date(argv[4]); p.jd_arrive = parse_date(argv[5]);

    // Parse --via pairs
    for (int i = 6; i < argc; ++i) {
        if (std::string(argv[i]) == "--via" && i + 2 < argc) {
            ++i;
            while (i + 1 < argc && argv[i][0] != '-') {
                p.flybys.push_back({argv[i], parse_date(argv[i+1])});
                i += 2;
            }
            --i;
        }
    }

    auto result = evaluate_multi_flyby(p, bodies);
    if (!result.valid) {
        std::cerr << "Multi-flyby evaluation failed (Lambert did not converge)\n";
        return 1;
    }

    std::cout << "Multi-Flyby Trajectory: " << p.departure_body;
    for (const auto& fb : p.flybys) std::cout << " -> " << fb.body_name;
    std::cout << " -> " << p.arrival_body << "\n\n";

    std::cout << std::left << std::setw(12) << "From" << std::setw(12) << "To"
              << std::right << std::setw(14) << "Depart" << std::setw(14) << "Arrive"
              << std::setw(10) << "TOF(d)" << std::setw(10) << "dv(km/s)" << "\n";
    std::cout << std::string(72, '-') << "\n";

    for (const auto& leg : result.legs) {
        int y, m, d; double h;
        std::cout << std::left << std::setw(12) << leg.from_body << std::setw(12) << leg.to_body;
        jd_to_calendar(leg.jd_depart, y, m, d, h);
        std::cout << std::right << "  " << y << "-" << std::setfill('0')
                  << std::setw(2) << m << "-" << std::setw(2) << d << std::setfill(' ');
        jd_to_calendar(leg.jd_arrive, y, m, d, h);
        std::cout << "  " << y << "-" << std::setfill('0')
                  << std::setw(2) << m << "-" << std::setw(2) << d << std::setfill(' ');
        std::cout << std::fixed << std::setw(10) << std::setprecision(0) << leg.tof_days
                  << std::setw(10) << std::setprecision(3) << leg.dv << "\n";
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nDeparture dv: " << result.dv_departure << " km/s\n";
    std::cout << "Arrival dv:   " << result.dv_arrival << " km/s\n";
    std::cout << "Total dv:     " << result.dv_total << " km/s\n";
    std::cout << "Total TOF:    " << std::setprecision(0) << result.total_tof_days << " days\n";
    return 0;
}

static int cmd_lagrange(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: solar lagrange <primary> <secondary>\n"
                  << "Example: solar lagrange Sun Earth\n";
        return 1;
    }
    auto sys = cr3bp::make_system(argv[2], argv[3]);
    auto lps = cr3bp::compute_lagrange_points(sys.mu);

    std::cout << "CR3BP System: " << sys.name << "\n";
    std::cout << "Mass parameter mu = " << std::scientific << std::setprecision(6) << sys.mu << "\n";
    std::cout << "L_star = " << std::setprecision(1) << sys.L_star << " km\n";
    std::cout << "T_star = " << sys.T_star << " s ("
              << std::fixed << std::setprecision(2) << sys.T_star / constants::DAY << " days)\n\n";

    struct { const char* name; cr3bp::CR3BPState s; } pts[] = {
        {"L1", lps.L1}, {"L2", lps.L2}, {"L3", lps.L3}, {"L4", lps.L4}, {"L5", lps.L5}
    };

    std::cout << std::left << std::setw(5) << "Point"
              << std::right << std::setw(16) << "x (norm)"
              << std::setw(16) << "y (norm)"
              << std::setw(16) << "x (km)"
              << std::setw(16) << "y (km)"
              << std::setw(12) << "Jacobi" << "\n";
    std::cout << std::string(81, '-') << "\n";

    for (const auto& p : pts) {
        double C = cr3bp::jacobi_constant(p.s, sys.mu);
        auto phys = cr3bp::to_physical(p.s, sys);
        std::cout << std::left << std::setw(5) << p.name << std::right
                  << std::fixed << std::setprecision(6)
                  << std::setw(16) << p.s.x
                  << std::setw(16) << p.s.y
                  << std::scientific << std::setprecision(4)
                  << std::setw(16) << phys.pos.x
                  << std::setw(16) << phys.pos.y
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << C << "\n";
    }

    // Distance from secondary body for L1 and L2
    double l1_dist = std::fabs(lps.L1.x - (1.0 - sys.mu)) * sys.L_star;
    double l2_dist = std::fabs(lps.L2.x - (1.0 - sys.mu)) * sys.L_star;
    std::cout << "\nL1 distance from " << argv[3] << ": " << std::scientific << std::setprecision(4) << l1_dist << " km\n";
    std::cout << "L2 distance from " << argv[3] << ": " << l2_dist << " km\n";

    return 0;
}

static int cmd_halo(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: solar halo <primary> <secondary> <L1|L2> <Az_km>\n"
                  << "Example: solar halo Sun Earth 1 110000\n";
        return 1;
    }

    auto sys = cr3bp::make_system(argv[2], argv[3]);
    int lpoint = std::atoi(argv[4]);
    double Az_km = std::atof(argv[5]);
    double Az_norm = Az_km / sys.L_star;

    cr3bp::HaloOrbitParams params;
    params.lpoint = lpoint;
    params.Az = Az_norm;
    params.family = 1; // Northern

    std::cerr << "Computing Halo orbit: " << sys.name << " L" << lpoint
              << ", Az=" << Az_km << " km (" << std::scientific << Az_norm << " norm)...\n";

    auto orbit = cr3bp::compute_halo_orbit(sys, params);

    if (!orbit.converged) {
        std::cerr << "Differential correction did not converge after "
                  << orbit.iterations << " iterations.\n";
        // Still print the Richardson IC
        auto ic = cr3bp::richardson_halo_ic(sys.mu, lpoint, Az_norm);
        std::cout << "Richardson analytical IC (not converged):\n";
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  x=" << ic.x << " y=" << ic.y << " z=" << ic.z << "\n";
        std::cout << "  xdot=" << ic.xdot << " ydot=" << ic.ydot << " zdot=" << ic.zdot << "\n";
        return 1;
    }

    double period_days = orbit.period * sys.T_star / constants::DAY;
    auto ic_phys = cr3bp::to_physical(orbit.initial_state, sys);

    std::cout << std::fixed;
    std::cout << "Halo Orbit: " << sys.name << " L" << lpoint << "\n";
    std::cout << "  Converged in " << orbit.iterations << " iterations\n\n";
    std::cout << "  Period: " << std::setprecision(6) << orbit.period << " (normalized), "
              << std::setprecision(2) << period_days << " days\n";
    std::cout << "  Jacobi constant: " << std::setprecision(6) << orbit.jacobi_constant << "\n";
    std::cout << "  Stability index: " << std::setprecision(4) << orbit.stability_index << "\n\n";

    std::cout << "  Initial state (normalized):\n";
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "    x=" << orbit.initial_state.x << "  z=" << orbit.initial_state.z << "\n";
    std::cout << "    ydot=" << orbit.initial_state.ydot << "\n\n";

    std::cout << "  Initial state (physical):\n";
    std::cout << "    x=" << ic_phys.pos.x << " km  z=" << ic_phys.pos.z << " km\n";
    std::cout << "    ydot=" << ic_phys.vel.y << " km/s\n";

    // Output trajectory if --trajectory flag
    bool output_traj = has_flag(argc, argv, "--trajectory");
    if (output_traj && !orbit.trajectory.empty()) {
        std::cout << "\n# trajectory: t  x  y  z  (normalized)\n";
        double dt_out = orbit.period / orbit.trajectory.size();
        for (size_t i = 0; i < orbit.trajectory.size(); ++i) {
            const auto& s = orbit.trajectory[i];
            std::cout << std::scientific << std::setprecision(6)
                      << i * dt_out << "\t" << s.x << "\t" << s.y << "\t" << s.z << "\n";
        }
    }

    return 0;
}

static int cmd_network(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: solar network <YYYY-MM-DD>\n";
        return 1;
    }
    double jd = parse_date(argv[2]);
    auto bodies = load_solar_system();

    SolarNetwork net;
    net.add_node({"Earth", NodeType::PlanetOrbit});
    net.add_node({"Mars", NodeType::PlanetOrbit});
    net.add_node({"Venus", NodeType::PlanetOrbit});
    net.add_node({"Jupiter", NodeType::PlanetOrbit});
    net.add_node({"Saturn", NodeType::PlanetOrbit});

    std::cerr << "Generating transfer network for " << argv[2] << "...\n";
    net.generate_lambert_edges(bodies, jd, 30, 1500, 30);
    net.print_summary();

    // Find and print shortest paths
    std::cout << "\nShortest paths (min dv):\n";
    const char* targets[] = {"Mars", "Jupiter", "Saturn"};
    for (const char* dest : targets) {
        auto path = net.shortest_path("Earth", dest, PathMetric::MinDv);
        if (path.valid) {
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "  Earth -> " << std::left << std::setw(10) << dest << std::right
                      << " dv=" << std::setw(7) << path.total_dv << " km/s"
                      << " TOF=" << std::setprecision(0) << std::setw(6) << path.total_tof_days << " d";
            if (path.nodes.size() > 2) {
                std::cout << "  via:";
                for (size_t i = 1; i < path.nodes.size() - 1; ++i)
                    std::cout << " " << path.nodes[i];
            }
            std::cout << "\n";
        }
    }

    return 0;
}

static int cmd_montecarlo(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: solar montecarlo <template> [--samples N]\n";
        return 1;
    }
    std::string tmpl = argv[2];
    MissionPlan plan;
    if (tmpl == "mars") plan = mars_direct_template();
    else if (tmpl == "gateway") plan = lunar_gateway_template();
    else if (tmpl == "grand-tour") plan = grand_tour_template();
    else { std::cerr << "Unknown template: " << tmpl << "\n"; return 1; }

    MCParams params;
    params.n_samples = 200;
    for (int i = 3; i < argc - 1; ++i) {
        if (std::string(argv[i]) == "--samples") params.n_samples = std::atoi(argv[i+1]);
    }

    std::cerr << "Running " << params.n_samples << " Monte Carlo samples for " << plan.name << "...\n";
    auto result = run_montecarlo(plan, params);
    print_mc_report(result, params);

    std::cerr << "\nRunning sensitivity analysis...\n";
    auto sens = run_sensitivity(plan);
    std::cout << "\nSensitivity (impact on success rate):\n";
    for (const auto& e : sens.entries) {
        std::cout << "  " << std::left << std::setw(22) << e.param_name
                  << std::right << std::fixed << std::setprecision(3)
                  << e.sensitivity * 100 << " pp\n";
    }

    return 0;
}

static int cmd_mission(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: solar mission <template>\n"
                  << "Templates: mars, gateway, grand-tour\n";
        return 1;
    }

    std::string tmpl = argv[2];
    MissionPlan plan;

    if (tmpl == "mars") plan = mars_direct_template();
    else if (tmpl == "gateway") plan = lunar_gateway_template();
    else if (tmpl == "grand-tour") plan = grand_tour_template();
    else {
        std::cerr << "Unknown template: " << tmpl << "\n";
        std::cerr << "Available: mars, gateway, grand-tour\n";
        return 1;
    }

    auto report = execute_mission(plan);
    print_mission_report(report);
    return report.success ? 0 : 1;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    // Parse global --de flag (must come before command)
    // Usage: solar --de <path> <command> [args...]
    // Also check SOLAR_DE_FILE environment variable
    static JPLEphemeris de_eph;
    std::string de_path;
    int cmd_start = 1; // index of the command argument

    if (argc >= 3 && std::string(argv[1]) == "--de") {
        de_path = argv[2];
        cmd_start = 3;
    } else {
        const char* env = std::getenv("SOLAR_DE_FILE");
        if (env) de_path = env;
    }

    if (!de_path.empty()) {
        if (de_eph.load_ascii(de_path)) {
            set_global_de_ephemeris(&de_eph);
            std::cerr << "DE ephemeris loaded: "
                      << de_eph.num_records() << " records, JD "
                      << std::fixed << std::setprecision(1)
                      << de_eph.jd_start() << " - " << de_eph.jd_end()
                      << ", EMRAT=" << std::setprecision(4) << de_eph.emrat() << "\n";
        } else {
            std::cerr << "Warning: failed to load DE file: " << de_path << "\n";
        }
    }

    if (cmd_start >= argc) {
        print_usage();
        return 1;
    }

    std::string cmd = argv[cmd_start];

    // Shift argv so subcommands see themselves at argv[1]
    int sub_argc = argc - cmd_start + 1;
    // Build shifted argv: argv[0] stays, then argv[cmd_start..]
    std::vector<char*> sub_argv;
    sub_argv.push_back(argv[0]);
    for (int i = cmd_start; i < argc; ++i) {
        sub_argv.push_back(argv[i]);
    }

    try {
        if (cmd == "bodies")     return cmd_bodies(sub_argc, sub_argv.data());
        if (cmd == "info")       return cmd_info(sub_argc, sub_argv.data());
        if (cmd == "moons")      return cmd_moons(sub_argc, sub_argv.data());
        if (cmd == "ephemeris")  return cmd_ephemeris(sub_argc, sub_argv.data());
        if (cmd == "orbits")     return cmd_orbits(sub_argc, sub_argv.data());
        if (cmd == "simulate")   return cmd_simulate(sub_argc, sub_argv.data());
        if (cmd == "transfer")   return cmd_transfer(sub_argc, sub_argv.data());
        if (cmd == "lambert")    return cmd_lambert(sub_argc, sub_argv.data());
        if (cmd == "energy")     return cmd_energy(sub_argc, sub_argv.data());
        if (cmd == "frame")     return cmd_frame(sub_argc, sub_argv.data());
        if (cmd == "time")      return cmd_time(sub_argc, sub_argv.data());
        if (cmd == "porkchop")  return cmd_porkchop(sub_argc, sub_argv.data());
        if (cmd == "launch-window") return cmd_launch_window(sub_argc, sub_argv.data());
        if (cmd == "flyby")     return cmd_flyby(sub_argc, sub_argv.data());
        if (cmd == "multi-flyby") return cmd_multi_flyby(sub_argc, sub_argv.data());
        if (cmd == "mission")  return cmd_mission(sub_argc, sub_argv.data());
        if (cmd == "montecarlo") return cmd_montecarlo(sub_argc, sub_argv.data());
        if (cmd == "network")  return cmd_network(sub_argc, sub_argv.data());
        if (cmd == "lagrange") return cmd_lagrange(sub_argc, sub_argv.data());
        if (cmd == "halo")     return cmd_halo(sub_argc, sub_argv.data());
        if (cmd == "blackhole") return run_blackhole_command(sub_argc, sub_argv.data());

        std::cerr << "Unknown command: " << cmd << "\n";
        print_usage();
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
