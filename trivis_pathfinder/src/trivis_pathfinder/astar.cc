/**
 * File:   astar.cc
 *
 * Date:   22.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_pathfinder/astar.h"

#include <iostream>

#include "boost/graph/astar_search.hpp"

using namespace trivis_pathfinder;

class TargetVisitor : public boost::default_astar_visitor {
public:
    struct TargetFound : public std::runtime_error {
        TargetFound() : std::runtime_error("target found") {}
    };
    explicit TargetVisitor(int target_graph_v) : _target_graph_v(target_graph_v) {}
    void examine_vertex(int v, const BoostGraph &) const {
        if (v == _target_graph_v) throw TargetFound();
    }
private:
    int _target_graph_v;
};

double trivis_pathfinder::ComputeCityPairShortestPathAstar(
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoint &target_city,
    bool cities_visible,
    const std::vector<int> &source_visible_reflex_vertices,
    const std::vector<int> &target_visible_reflex_vertices,
    const trivis::geom::FPoints &reflex_vertices_points,
    BoostGraph &reflex_graph,
    std::vector<int> &id_path_no_cities,
    bool revert_graph
) {
    if (cities_visible) {
        id_path_no_cities.clear();
        return source_city.DistanceTo(target_city);
    }
    if (source_city == target_city) {
        id_path_no_cities.clear();
        return 0.0;
    }
    int n_reflex_vertices = static_cast<int>(reflex_vertices_points.size());
    int n_graph_orig = static_cast<int>(boost::num_vertices(reflex_graph));

    // Add vertex s and its edges.
    int source_graph_v = n_graph_orig;
    boost::add_vertex(reflex_graph);
    for (int reflex_vertex_id: source_visible_reflex_vertices) {
        const auto &reflex_vertex_point = reflex_vertices_points[reflex_vertex_id];
        double dist = source_city.DistanceTo(reflex_vertex_point);
        boost::add_edge(source_graph_v, reflex_vertex_id, dist, reflex_graph);
    }

    // Add vertex g and its edges.
    boost::add_vertex(reflex_graph);
    int target_graph_v = n_graph_orig + 1;
    for (int reflex_vertex_id: target_visible_reflex_vertices) {
        const auto &reflex_vertex_point = reflex_vertices_points[reflex_vertex_id];
        double dist = target_city.DistanceTo(reflex_vertex_point);
        boost::add_edge(reflex_vertex_id, target_graph_v, dist, reflex_graph);
    }

    // Create the heuristic.
    auto heuristic = [source_graph_v, target_graph_v, &source_city, &target_city, n_reflex_vertices, &reflex_vertices_points](int graph_v) {
        if (graph_v == source_graph_v) return source_city.DistanceTo(target_city);
        if (graph_v == target_graph_v) return 0.0;
        if (graph_v >= n_reflex_vertices) return std::numeric_limits<double>::max();
        return reflex_vertices_points[graph_v].DistanceTo(target_city);
    };

    // Run AStar.
    std::vector<int> p(n_graph_orig + 2);
    std::vector<double> d(n_graph_orig + 2);
    try {
        auto param = boost::predecessor_map(&p[0]).distance_map(&d[0]).visitor(TargetVisitor(target_graph_v));
        boost::astar_search(reflex_graph, source_graph_v, heuristic, param);
    } catch (const TargetVisitor::TargetFound &e) {
        // Catching this exception means Astar found the target.
        // This is desired behavior: do nothing.
    }

    // Reconstruct the path.
    id_path_no_cities.clear();
    int graph_v = target_graph_v;
    do {
        graph_v = p[graph_v];
        id_path_no_cities.push_back(graph_v);
    } while (graph_v != source_graph_v);
    id_path_no_cities.pop_back(); // remove source

    // Reverse the path to get the correct order.
    std::reverse(id_path_no_cities.begin(), id_path_no_cities.end());

    if (revert_graph) {

        // Return the graph to the original state.
        boost::clear_vertex(source_graph_v, reflex_graph);
        boost::clear_vertex(target_graph_v, reflex_graph);
        boost::remove_vertex(source_graph_v, reflex_graph);
        boost::remove_vertex(target_graph_v, reflex_graph);

    }

    return d[target_graph_v];
}
