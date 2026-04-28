#include "solar/network.h"
#include "solar/ephemeris.h"
#include "solar/transfer.h"
#include "solar/constants.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>
#include <iomanip>

namespace solar {

void SolarNetwork::add_node(const NetworkNode& node) {
    nodes_.push_back(node);
}

int SolarNetwork::find_node(const std::string& name) const {
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].name == name) return static_cast<int>(i);
    }
    return -1;
}

void SolarNetwork::generate_lambert_edges(
    const std::vector<Body>& bodies,
    double jd_epoch,
    double tof_min_days, double tof_max_days,
    int tof_steps)
{
    double dt = (tof_max_days - tof_min_days) / std::max(1, tof_steps - 1);

    // For each pair of planet-orbit nodes, find best Lambert transfer
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].type != NodeType::PlanetOrbit) continue;
        const Body* b_from = find_body(bodies, nodes_[i].name);
        if (!b_from) continue;

        for (size_t j = 0; j < nodes_.size(); ++j) {
            if (i == j) continue;
            if (nodes_[j].type != NodeType::PlanetOrbit) continue;
            const Body* b_to = find_body(bodies, nodes_[j].name);
            if (!b_to) continue;

            NetworkEdge best;
            best.from = static_cast<int>(i);
            best.to = static_cast<int>(j);
            best.dv_total = std::numeric_limits<double>::max();
            best.valid = false;

            State s_from = compute_position(*b_from, jd_epoch);

            for (int k = 0; k < tof_steps; ++k) {
                double tof_days = tof_min_days + k * dt;
                double jd_arrive = jd_epoch + tof_days;
                double tof_sec = tof_days * constants::DAY;

                State s_to = compute_position(*b_to, jd_arrive);

                auto lam = solve_lambert(s_from.pos, s_to.pos, tof_sec,
                                         constants::MU_SUN, true);
                if (!lam.converged) continue;

                Vec3 dv1_vec = lam.v1 - s_from.vel;
                Vec3 dv2_vec = s_to.vel - lam.v2;
                double dv1 = dv1_vec.norm();
                double dv2 = dv2_vec.norm();
                double dv_total = dv1 + dv2;

                if (dv_total < best.dv_total) {
                    best.dv_depart = dv1;
                    best.dv_arrive = dv2;
                    best.dv_total = dv_total;
                    best.tof_days = tof_days;
                    best.jd_depart = jd_epoch;
                    best.valid = true;
                }
            }

            if (best.valid) {
                edges_.push_back(best);
            }
        }
    }
}

NetworkEdge SolarNetwork::get_best_edge(const std::string& from,
                                         const std::string& to) const {
    int fi = find_node(from);
    int ti = find_node(to);
    NetworkEdge best;
    best.valid = false;
    best.dv_total = std::numeric_limits<double>::max();

    for (const auto& e : edges_) {
        if (e.from == fi && e.to == ti && e.dv_total < best.dv_total) {
            best = e;
        }
    }
    return best;
}

PathResult SolarNetwork::shortest_path(const std::string& from,
                                        const std::string& to,
                                        PathMetric metric) const {
    PathResult result;
    result.valid = false;

    int src = find_node(from);
    int dst = find_node(to);
    if (src < 0 || dst < 0) return result;

    int n = static_cast<int>(nodes_.size());
    std::vector<double> dist(n, std::numeric_limits<double>::max());
    std::vector<int> prev(n, -1);
    std::vector<int> prev_edge(n, -1);

    // Build adjacency list
    std::vector<std::vector<int>> adj(n);
    for (size_t i = 0; i < edges_.size(); ++i) {
        adj[edges_[i].from].push_back(static_cast<int>(i));
    }

    // Dijkstra
    using PQEntry = std::pair<double, int>; // (cost, node)
    std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<>> pq;

    dist[src] = 0;
    pq.push({0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;
        if (u == dst) break;

        for (int ei : adj[u]) {
            const auto& e = edges_[ei];
            double w = (metric == PathMetric::MinDv) ? e.dv_total : e.tof_days;
            double nd = dist[u] + w;
            if (nd < dist[e.to]) {
                dist[e.to] = nd;
                prev[e.to] = u;
                prev_edge[e.to] = ei;
                pq.push({nd, e.to});
            }
        }
    }

    if (dist[dst] >= std::numeric_limits<double>::max()) return result;

    // Reconstruct path
    std::vector<int> path_nodes;
    for (int v = dst; v != -1; v = prev[v]) {
        path_nodes.push_back(v);
    }
    std::reverse(path_nodes.begin(), path_nodes.end());

    result.valid = true;
    result.total_dv = 0;
    result.total_tof_days = 0;

    for (int v : path_nodes) {
        result.nodes.push_back(nodes_[v].name);
    }

    for (size_t i = 1; i < path_nodes.size(); ++i) {
        int ei = prev_edge[path_nodes[i]];
        if (ei >= 0) {
            result.edges.push_back(edges_[ei]);
            result.total_dv += edges_[ei].dv_total;
            result.total_tof_days += edges_[ei].tof_days;
        }
    }

    return result;
}

void SolarNetwork::print_summary() const {
    std::cout << "Solar Network: " << nodes_.size() << " nodes, "
              << edges_.size() << " edges\n";
    std::cout << "\nNodes:\n";
    for (const auto& n : nodes_) {
        const char* type = "?";
        switch (n.type) {
            case NodeType::PlanetOrbit: type = "orbit"; break;
            case NodeType::LagrangePoint: type = "L-point"; break;
            case NodeType::Depot: type = "depot"; break;
            case NodeType::Asteroid: type = "asteroid"; break;
        }
        std::cout << "  " << std::left << std::setw(12) << n.name << type << "\n";
    }
    std::cout << "\nEdges (best transfer per pair):\n";
    std::cout << std::left << std::setw(12) << "From" << std::setw(12) << "To"
              << std::right << std::setw(10) << "dv(km/s)" << std::setw(10) << "TOF(d)" << "\n";
    std::cout << std::string(44, '-') << "\n";
    for (const auto& e : edges_) {
        std::cout << std::left << std::setw(12) << nodes_[e.from].name
                  << std::setw(12) << nodes_[e.to].name
                  << std::right << std::fixed
                  << std::setw(10) << std::setprecision(3) << e.dv_total
                  << std::setw(10) << std::setprecision(0) << e.tof_days << "\n";
    }
}

} // namespace solar
