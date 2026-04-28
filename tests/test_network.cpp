// TEST-FIRST: Solar system transportation network
//
// Acceptance criteria:
// 1. SolarNetwork holds nodes and edges with dv/time/risk weights
// 2. Can auto-generate Lambert edges between planet nodes for a given epoch
// 3. Dijkstra finds minimum-dv path
// 4. Can compare direct route vs via-staging-node route
// 5. Outputs real numbers that make physical sense

#include "solar/network.h"
#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace solar;

static int pass = 0, fail_count = 0;

static void check(const char* name, bool condition) {
    if (condition) { std::cout << "  PASS: " << name << "\n"; pass++; }
    else { std::cout << "  FAIL: " << name << "\n"; fail_count++; }
}

int main() {
    std::cout << "=== Solar Network Acceptance Tests ===\n\n";

    // --- Test 1: Build a network with nodes ---
    std::cout << "Test 1: Create network with nodes\n";
    {
        SolarNetwork net;
        net.add_node({"Earth", NodeType::PlanetOrbit});
        net.add_node({"Mars", NodeType::PlanetOrbit});
        net.add_node({"Venus", NodeType::PlanetOrbit});
        net.add_node({"Jupiter", NodeType::PlanetOrbit});

        check("has 4 nodes", net.num_nodes() == 4);
        check("can find Earth", net.find_node("Earth") >= 0);
        check("unknown returns -1", net.find_node("Pluto") == -1);
    }

    // --- Test 2: Auto-generate Lambert edges for 2026 window ---
    std::cout << "\nTest 2: Generate Lambert edges from ephemeris\n";
    {
        SolarNetwork net;
        auto bodies = load_solar_system();

        // Add inner planet nodes
        net.add_node({"Earth", NodeType::PlanetOrbit});
        net.add_node({"Mars", NodeType::PlanetOrbit});
        net.add_node({"Venus", NodeType::PlanetOrbit});

        // Generate edges for a specific epoch window
        double jd = calendar_to_jd(2026, 10, 1);
        net.generate_lambert_edges(bodies, jd, 100, 400); // TOF range: 100-400 days

        check("edges generated", net.num_edges() > 0);

        // Earth-Mars edge should exist with dv ~5-7 km/s
        auto em = net.get_best_edge("Earth", "Mars");
        check("Earth-Mars edge exists", em.valid);
        check("Earth-Mars dv reasonable", em.dv_total > 4.0 && em.dv_total < 10.0);
        check("Earth-Mars tof reasonable", em.tof_days > 100 && em.tof_days < 400);

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Earth-Mars: dv=" << em.dv_total << " km/s, TOF="
                  << std::setprecision(0) << em.tof_days << " days\n";
    }

    // --- Test 3: Dijkstra shortest path ---
    std::cout << "\nTest 3: Find minimum-dv path\n";
    {
        SolarNetwork net;
        auto bodies = load_solar_system();

        net.add_node({"Earth", NodeType::PlanetOrbit});
        net.add_node({"Mars", NodeType::PlanetOrbit});
        net.add_node({"Venus", NodeType::PlanetOrbit});
        net.add_node({"Jupiter", NodeType::PlanetOrbit});

        double jd = calendar_to_jd(2026, 10, 1);
        net.generate_lambert_edges(bodies, jd, 80, 500);

        auto path = net.shortest_path("Earth", "Mars", PathMetric::MinDv);

        check("path found", path.valid);
        check("path starts at Earth", !path.nodes.empty() && path.nodes.front() == "Earth");
        check("path ends at Mars", !path.nodes.empty() && path.nodes.back() == "Mars");
        check("total dv > 0", path.total_dv > 0);

        std::cout << "  Path: ";
        for (size_t i = 0; i < path.nodes.size(); ++i) {
            if (i > 0) std::cout << " -> ";
            std::cout << path.nodes[i];
        }
        std::cout << "\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Total dv: " << path.total_dv << " km/s, TOF: "
                  << std::setprecision(0) << path.total_tof_days << " days\n";
    }

    // --- Test 4: Compare direct vs via-Venus ---
    std::cout << "\nTest 4: Route comparison (direct vs staging)\n";
    {
        SolarNetwork net;
        auto bodies = load_solar_system();

        net.add_node({"Earth", NodeType::PlanetOrbit});
        net.add_node({"Mars", NodeType::PlanetOrbit});
        net.add_node({"Venus", NodeType::PlanetOrbit});

        double jd = calendar_to_jd(2026, 10, 1);
        net.generate_lambert_edges(bodies, jd, 80, 400);

        auto direct = net.shortest_path("Earth", "Mars", PathMetric::MinDv);
        // Force via Venus by checking if Venus route exists
        auto via_venus = net.shortest_path("Earth", "Venus", PathMetric::MinDv);
        auto venus_mars = net.shortest_path("Venus", "Mars", PathMetric::MinDv);

        if (direct.valid) {
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "  Direct Earth-Mars: " << direct.total_dv << " km/s\n";
        }
        if (via_venus.valid && venus_mars.valid) {
            double via_dv = via_venus.total_dv + venus_mars.total_dv;
            std::cout << "  Via Venus: " << via_dv << " km/s (E-V: "
                      << via_venus.total_dv << " + V-M: " << venus_mars.total_dv << ")\n";
        }

        check("direct route found", direct.valid);
        check("comparison possible", via_venus.valid && venus_mars.valid);
    }

    std::cout << "\n=== Results: " << pass << " passed, " << fail_count << " failed ===\n";
    return fail_count > 0 ? 1 : 0;
}
