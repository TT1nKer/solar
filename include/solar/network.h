#pragma once

#include "solar/body.h"
#include <vector>
#include <string>

namespace solar {

enum class NodeType { PlanetOrbit, LagrangePoint, Depot, Asteroid };
enum class PathMetric { MinDv, MinTime };

struct NetworkNode {
    std::string name;
    NodeType type;
};

struct NetworkEdge {
    int from, to;          // node indices
    double dv_total;       // km/s (departure + arrival)
    double dv_depart;
    double dv_arrive;
    double tof_days;
    double jd_depart;
    bool valid = false;
};

struct PathResult {
    std::vector<std::string> nodes;   // node names in order
    std::vector<NetworkEdge> edges;   // edges along the path
    double total_dv;
    double total_tof_days;
    bool valid = false;
};

class SolarNetwork {
public:
    void add_node(const NetworkNode& node);
    int find_node(const std::string& name) const;
    int num_nodes() const { return static_cast<int>(nodes_.size()); }
    int num_edges() const { return static_cast<int>(edges_.size()); }

    // Generate Lambert transfer edges between all planet-orbit nodes
    // tof_min/tof_max: range of transfer times to scan (days)
    // Picks the best (min dv) transfer for each pair
    void generate_lambert_edges(const std::vector<Body>& bodies,
                                double jd_epoch,
                                double tof_min_days, double tof_max_days,
                                int tof_steps = 20);

    // Get the best direct edge between two nodes
    NetworkEdge get_best_edge(const std::string& from, const std::string& to) const;

    // Dijkstra shortest path
    PathResult shortest_path(const std::string& from, const std::string& to,
                             PathMetric metric = PathMetric::MinDv) const;

    // Print network summary
    void print_summary() const;

private:
    std::vector<NetworkNode> nodes_;
    std::vector<NetworkEdge> edges_;
};

} // namespace solar
