/**
 * File:   dijkstra.cc
 *
 * Date:   21.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_pathfinder/algorithms/dijkstra.h"

#include <iostream>

#include "boost/graph/dijkstra_shortest_paths.hpp"

#include "trivis/trivis.h"

using namespace trivis_pathfinder;
using namespace trivis_pathfinder::algorithms;

class EarlyExitVisitor : public boost::default_bfs_visitor {
public:
    struct Finished : public std::runtime_error {
        Finished() : std::runtime_error("finished") {}
    };
    EarlyExitVisitor(int n_reflex_vertices, int city_cnt_max) : _n_reflex_vertices(n_reflex_vertices), _city_cnt_max(city_cnt_max) {}
    void initialize_vertex(int, const BoostGraph &) const {}
    void discover_vertex(int, const BoostGraph &) const {}
    void examine_vertex(int graph_v, const BoostGraph &) {
        if (graph_v >= _n_reflex_vertices) {
            ++_city_cnt; // count cities
            if (_city_cnt >= _city_cnt_max) {
                // finish when all cities are visited
                throw Finished();
            }
        }
    }
    void examine_edge(const BoostGraphEdgeDescriptor &, const BoostGraph &) const {}
    void edge_relaxed(const BoostGraphEdgeDescriptor &, const BoostGraph &) const {}
    void edge_not_relaxed(const BoostGraphEdgeDescriptor &, const BoostGraph &) const {}
    void finish_vertex(int, const BoostGraph &) const {}
protected:
    const int _n_reflex_vertices;
    const int _city_cnt_max;
    int _city_cnt = 0;
};

void algorithms::ComputeShortestPathsDijkstraAllPairsReflex(
    const BoostGraph &reflex_graph,
    std::vector<std::vector<double>> &paths_lengths,
    std::vector<std::vector<int>> &paths_predecessor_map
) {
    int n_reflex_vertices = static_cast<int>(boost::num_vertices(reflex_graph));
    paths_lengths = std::vector<std::vector<double>>(n_reflex_vertices, std::vector<double>(n_reflex_vertices, -1.0));
    paths_predecessor_map = std::vector<std::vector<int>>(n_reflex_vertices, std::vector<int>(n_reflex_vertices, -1));
    for (int source_id = 0; source_id < n_reflex_vertices; ++source_id) {
        // Find the shortest paths.
        auto &d = paths_lengths[source_id];
        auto &p = paths_predecessor_map[source_id];
        boost::dijkstra_shortest_paths(reflex_graph, source_id, boost::predecessor_map(&p[0]).distance_map(&d[0]));
    }
}

std::vector<int> algorithms::GetShortestIDPathAllPairsReflex(
    int source_id,
    int target_id,
    const std::vector<int> &source_predecessors
) {
    if (source_id == target_id) {
        return std::vector<int>{};
    }
    std::vector<int> id_path;
    int source_graph_v = source_id;
    int curr_graph_v = target_id;
    id_path.push_back(curr_graph_v);
    do {
        curr_graph_v = source_predecessors[curr_graph_v];
        id_path.push_back(curr_graph_v);
    } while (curr_graph_v != source_graph_v);
    // Now path starts in target and ends in source.
    // Reverse the path to get the correct order.
    std::reverse(id_path.begin(), id_path.end());
    return id_path;
}

double algorithms::ComputeShortestPathDijkstraCities(
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoint &target_city,
    bool cities_visible,
    const std::vector<int> &source_visible_reflex_vertices,
    const std::vector<int> &target_visible_reflex_vertices,
    const trivis::geom::FPoints &reflex_vertices_points,
    const std::vector<std::vector<double>> &reflex_vertices_paths_lengths,
    const std::vector<std::vector<int>> &reflex_vertices_paths_predecessor_map,
    std::vector<int> &id_path_no_endpoints
) {
    if (cities_visible) {
        id_path_no_endpoints.clear();
        return source_city.DistanceTo(target_city);
    }
    double min_length = std::numeric_limits<double>::max();
    int best_reflex_source_id = -1;
    int best_reflex_target_id = -1;
    for (int reflex_source_id: source_visible_reflex_vertices) {
        const auto &reflex_source_point = reflex_vertices_points[reflex_source_id];
        const auto &reflex_source_lengths = reflex_vertices_paths_lengths[reflex_source_id];
        double dist_i = source_city.DistanceTo(reflex_source_point);
        if (dist_i > min_length) {
            continue;
        }
        for (int reflex_target_id: target_visible_reflex_vertices) {
            const auto &reflex_target_point = reflex_vertices_points[reflex_target_id];
            double length = dist_i + reflex_source_lengths[reflex_target_id];
            if (length > min_length) {
                continue;
            }
            length += target_city.DistanceTo(reflex_target_point);
            if (length < min_length) {
                min_length = length;
                best_reflex_source_id = reflex_source_id;
                best_reflex_target_id = reflex_target_id;
            }
        }
    }
    if (best_reflex_source_id == -1 || best_reflex_target_id == -1) {
        id_path_no_endpoints.clear();
        return -1.0;
    }
    if (best_reflex_source_id == best_reflex_target_id) {
        id_path_no_endpoints = {best_reflex_source_id};
    } else {
        id_path_no_endpoints = GetShortestIDPathAllPairsReflex(best_reflex_source_id, best_reflex_target_id, reflex_vertices_paths_predecessor_map[best_reflex_source_id]);
    }
    return min_length;
}

trivis::geom::FPoints algorithms::ConvertIDPathToPointPathCities(
    const std::vector<int> &id_path_no_endpoints,
    const trivis::geom::FPoints &reflex_vertices_points,
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoint &target_city
) {
    auto point_path = trivis::utils::Concatenate({source_city}, trivis::utils::Select(reflex_vertices_points, id_path_no_endpoints));
    point_path.push_back(target_city);
    return point_path;
}

void algorithms::ComputeShortestPathsDijkstraAllPairsCities(
    BoostGraph reflex_graph,
    const trivis::geom::FPoints &reflex_vertices_points,
    const trivis::geom::FPoints &cities,
    const std::vector<std::vector<int>> &vis_graph_city_city,
    const std::vector<std::vector<int>> &vis_graph_city_reflex,
    std::vector<std::vector<double>> &paths_lengths,
    std::vector<std::vector<int>> &paths_predecessor_map
) {

    int n_reflex_vertices = static_cast<int>(reflex_vertices_points.size());
    int n_cities = static_cast<int>(cities.size());
    int n_graph = n_reflex_vertices + n_cities;

    paths_lengths = std::vector<std::vector<double>>(n_cities, std::vector<double>(n_graph, -1.0));
    paths_predecessor_map = std::vector<std::vector<int>>(n_cities, std::vector<int>(n_graph, -1));

    // Resize graph vertices.
    reflex_graph.m_vertices.resize(n_graph);
    for (int source_city_id = 0; source_city_id < n_cities; ++source_city_id) {
        // Add city vertex.
        int source_graph_v = n_reflex_vertices + source_city_id;
        reflex_graph.added_vertex(source_graph_v);
        // Add city-vertex edges (gradually).
        for (int reflex_vertex_v: vis_graph_city_reflex[source_city_id]) {
            double distance = cities[source_city_id].DistanceTo(reflex_vertices_points[reflex_vertex_v]);
            boost::add_edge(source_graph_v, reflex_vertex_v, distance, reflex_graph);
        }
        // Add city-city edges (gradually).
        for (int target_city_id: vis_graph_city_city[source_city_id]) {
            if (target_city_id >= source_city_id) {
                continue;
            }
            int target_graph_v = n_reflex_vertices + target_city_id;
            double distance = cities[source_city_id].DistanceTo(cities[target_city_id]);
            boost::add_edge(source_graph_v, target_graph_v, distance, reflex_graph);
        }
        if (source_city_id == 0) {
            paths_lengths[0][0] = 0.0;
            continue;
        }
        // Find the shortest paths.
        auto &d = paths_lengths[source_city_id];
        auto &p = paths_predecessor_map[source_city_id];
        try {
            auto param = boost::predecessor_map(&p[0]).distance_map(&d[0]).visitor(EarlyExitVisitor(n_reflex_vertices, source_city_id + 1));
            boost::dijkstra_shortest_paths(reflex_graph, source_graph_v, param);
            // Note: predecessor_map and distance_map will be valid only for s_id (start vertex) and g_id (destination vertex) where s_id > g_id (for efficiency).
            // To get the shortest path/dist for s_id < g_id we must use the shortest paths' symmetry.
        } catch (const EarlyExitVisitor::Finished &e) {
            // Catching this exception means Dijkstra finished prematurely.
            // This is desired behavior: do nothing.
        }
    }
}

double algorithms::GetShortestPathLengthAllPairsCities(
    int source_city_id,
    int target_city_id,
    int n_reflex_vertices,
    const std::vector<std::vector<double>> &paths_lengths
) {
    if (source_city_id == target_city_id) {
        return 0.0;
    }
    if (source_city_id < target_city_id) {
        return GetShortestPathLengthAllPairsCities(target_city_id, source_city_id, n_reflex_vertices, paths_lengths);
    }
    return paths_lengths[source_city_id][n_reflex_vertices + target_city_id];
}

std::vector<int> algorithms::GetShortestIDPathAllPairsCities(
    int source_city_id,
    int target_city_id,
    int n_reflex_vertices,
    const std::vector<std::vector<int>> &paths_predecessor_map,
    bool reversed
) {
    if (source_city_id == target_city_id) {
        return std::vector<int>{};
    }
    if (source_city_id < target_city_id) {
        return GetShortestIDPathAllPairsCities(target_city_id, source_city_id, n_reflex_vertices, paths_predecessor_map, !reversed);
    }
    std::vector<int> id_path;
    const auto &predecessors = paths_predecessor_map[source_city_id];
    int source_graph_v = n_reflex_vertices + source_city_id;
    int target_graph_v = n_reflex_vertices + target_city_id;
    id_path.push_back(target_graph_v);
    int curr_graph_v = target_graph_v;
    do {
        curr_graph_v = predecessors[curr_graph_v];
        id_path.push_back(curr_graph_v);
    } while (curr_graph_v != source_graph_v);
    // Now path starts from target and ends in source.
    // Reverse the path to get the correct order, unless reversed flag is true, then the reversed order is preserved.
    if (!reversed) {
        std::reverse(id_path.begin(), id_path.end());
    }
    return id_path;
}

trivis::geom::FPoints algorithms::ConvertIDPathToPointPathAllPairsCities(
    const std::vector<int> &id_path,
    const trivis::geom::FPoints &reflex_vertices_points,
    const trivis::geom::FPoints &cities
) {
    int n_reflex_vertices = static_cast<int>(reflex_vertices_points.size());
    int n_cities = static_cast<int>(cities.size());;
    int n_graph = n_reflex_vertices + n_cities;
    int n_path = static_cast<int>(id_path.size());
    trivis::geom::FPoints point_path;
    point_path.reserve(n_path);
    for (int graph_v: id_path) {
        if (graph_v < n_reflex_vertices) {
            point_path.push_back(reflex_vertices_points[graph_v]);
        } else if (graph_v < n_graph) {
            point_path.push_back(cities[graph_v - n_reflex_vertices]);
        }
    }
    return point_path;
}
