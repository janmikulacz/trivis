/**
 * File:   visibility_planner.cc
 *
 * Date:   24.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_pathfinder/visibility_planner.h"

#include <boost/heap/fibonacci_heap.hpp>

using namespace trivis_pathfinder;

void VisibilityPlanner::Init(
    const std::vector<bool> &is_node_convex,
    const std::vector<int> &mesh_id_to_reflex_id,
    const trivis::geom::FPoints &reflex_points,
    const std::vector<std::vector<int>> &vis_graph_reflex_reflex
) {
    _is_node_convex = is_node_convex;
    _mesh_id_to_reflex_id = mesh_id_to_reflex_id;
    _n_reflex = static_cast<int>(reflex_points.size());
    _source_v = _n_reflex;
    _target_v = _n_reflex + 1;
    _n_graph = _n_reflex + 2;
    _graph_neigh = std::vector<std::vector<int>>(_n_graph);
    _graph_costs = std::vector<std::vector<double>>(_n_graph);
    _graph_points = trivis::geom::FPoints(reflex_points);
    _d = std::vector<double>(_n_graph, std::numeric_limits<double>::max());
    _c = std::vector<Color>(_n_graph, Color::kWhite);
    _p = std::vector<int>(_n_graph, 0);
    for (int rfx_ver_id = 0; rfx_ver_id < _n_reflex; ++rfx_ver_id) {
        for (int visible_rfx_ver_id: vis_graph_reflex_reflex[rfx_ver_id]) {
            AddEdge(rfx_ver_id, visible_rfx_ver_id);
        }
    }
    _graph_points.resize(_n_graph);
}

double VisibilityPlanner::FindShortestPath(
    int source_reflex_vertex_id,
    int target_reflex_vertex_id,
    std::vector<int> &id_path
) {
    return FindShortestPath(id_path, source_reflex_vertex_id, target_reflex_vertex_id);
}

double VisibilityPlanner::FindShortestPath(
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoint &target_city,
    const trivis::Trivis &vis,
    std::vector<int> &id_path_no_endpoints
) {

    auto source_city_pl = vis.LocatePoint(source_city);
    if (!source_city_pl.has_value()) {
        //LOGF_ERR("[VisibilityPlanner::FindShortestPath] Source out of map!");
        id_path_no_endpoints.clear();
        return -1.0;
    }

    if (vis.IsVisible(source_city, source_city_pl.value(), target_city)) {
        // The target is directly visible from the source.
        id_path_no_endpoints.clear();
        return source_city.DistanceTo(target_city);
    }

    // Compute vertices visible from the source
    auto visible_vertices_source = vis.VisibleVertices(source_city, source_city_pl.value(), std::nullopt, _is_node_convex);

    auto target_city_pl = vis.LocatePoint(target_city);
    if (!target_city_pl.has_value()) {
        //LOGF_ERR("[VisibilityPlanner::FindShortestPath] Target out of map!");
        id_path_no_endpoints.clear();
        return -1.0;
    }

    auto visible_vertices_target = vis.VisibleVertices(target_city, target_city_pl.value(), std::nullopt, _is_node_convex);

    // Add the source to the graph.
    _graph_points[_source_v] = source_city;
    for (int vertex_id: visible_vertices_source) {
        AddEdge(_source_v, _mesh_id_to_reflex_id[vertex_id]);
    }

    // Add the target to the graph.
    _graph_points[_target_v] = target_city;
    for (int vertex_id: visible_vertices_target) {
        AddEdge(_mesh_id_to_reflex_id[vertex_id], _target_v);
    }

    double path_length = FindShortestPath(id_path_no_endpoints);

    // Cleaning after source and target.
    _graph_neigh[_source_v].clear();
    _graph_costs[_source_v].clear();
    for (auto vertex_id: visible_vertices_target) {
        _graph_neigh[_mesh_id_to_reflex_id[vertex_id]].pop_back();
        _graph_costs[_mesh_id_to_reflex_id[vertex_id]].pop_back();
    }

    return path_length;
}

double VisibilityPlanner::FindShortestPath(
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoint &target_city,
    bool target_visible,
    const std::vector<int> &source_visible_reflex_vertices,
    const std::vector<int> &target_visible_reflex_vertices,
    std::vector<int> &id_path_no_endpoints
) {

    if (target_visible) {
        // Target is directly visible from the source.
        id_path_no_endpoints.clear();
        return source_city.DistanceTo(target_city);
    }

    // Add the source to the graph.
    _graph_points[_source_v] = source_city;
    for (int reflex_vertex_id: source_visible_reflex_vertices) {
        AddEdge(_source_v, reflex_vertex_id);
    }

    // Add the target to the graph.
    _graph_points[_target_v] = target_city;
    for (int reflex_vertex_id: target_visible_reflex_vertices) {
        AddEdge(reflex_vertex_id, _target_v);
    }

    double path_length = FindShortestPath(id_path_no_endpoints);

    // Cleaning after source and target.
    _graph_neigh[_source_v].clear();
    _graph_costs[_source_v].clear();
    for (auto reflex_vertex_id: target_visible_reflex_vertices) {
        _graph_neigh[reflex_vertex_id].pop_back();
        _graph_costs[reflex_vertex_id].pop_back();
    }

    return path_length;
}

void VisibilityPlanner::FindShortestPaths(
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoints &target_cities,
    const trivis::Trivis &vis,
    std::vector<std::vector<int>> &id_paths_no_endpoints,
    std::vector<double> &paths_lengths
) {

    auto source_city_pl = vis.LocatePoint(source_city);
    if (!source_city_pl.has_value()) {
        //LOGF_ERR("[VisibilityPlanner::FindShortestPath] Source out of map!");
        id_paths_no_endpoints.clear();
        paths_lengths.clear();
        return;
    }

    // Compute vertices visible from the source
    auto visible_vertices_source = vis.VisibleVertices(source_city, source_city_pl.value(), std::nullopt, _is_node_convex);

    // Add the source to the graph.
    _graph_points[_source_v] = source_city;
    for (auto vertex_id: visible_vertices_source) {
        AddEdge(_source_v, _mesh_id_to_reflex_id[vertex_id]);
    }

    int n_targets = static_cast<int>(target_cities.size());
    id_paths_no_endpoints = std::vector<std::vector<int>>(n_targets);
    paths_lengths = std::vector<double>(n_targets, -1.0);

    for (int target_city_id = 0; target_city_id < n_targets; ++target_city_id) {
        const auto &target_city = target_cities[target_city_id];

        if (vis.IsVisible(source_city, source_city_pl.value(), target_city)) {
            // The target is directly visible from the source.
            paths_lengths[target_city_id] = source_city.DistanceTo(target_city);
            continue;
        }

        auto target_city_pl = vis.LocatePoint(target_city);
        if (!target_city_pl.has_value()) {
            //LOGF_WRN("[VisibilityPlanner::FindShortestPath] Target " << target_city_id << " out of map!");
            continue;
        }

        auto visible_vertices_target = vis.VisibleVertices(target_city, target_city_pl.value(), std::nullopt, _is_node_convex);

        // Add the target to the graph.
        _graph_points[_target_v] = target_city;
        for (auto vertex_id: visible_vertices_target) {
            AddEdge(_mesh_id_to_reflex_id[vertex_id], _target_v);
        }

        paths_lengths[target_city_id] = FindShortestPath(id_paths_no_endpoints[target_city_id]);

        // Cleaning after target.
        for (auto vertex_id: visible_vertices_target) {
            _graph_neigh[_mesh_id_to_reflex_id[vertex_id]].pop_back();
            _graph_costs[_mesh_id_to_reflex_id[vertex_id]].pop_back();
        }
    }
    // Cleaning after source.
    _graph_neigh[_source_v].clear();
    _graph_costs[_source_v].clear();
}

void VisibilityPlanner::FindShortestPaths(
    const trivis::geom::FPoint &source_city,
    const trivis::geom::FPoints &target_cities,
    std::vector<bool> &targets_visible,
    const std::vector<int> &source_visible_reflex_vertices,
    const std::vector<std::vector<int>> &targets_visible_reflex_vertices,
    std::vector<std::vector<int>> &id_paths_no_endpoints,
    std::vector<double> &paths_lengths
) {

    // Add the source to the graph.
    _graph_points[_source_v] = source_city;
    for (auto reflex_vertex_id: source_visible_reflex_vertices) {
        AddEdge(_source_v, reflex_vertex_id);
    }

    int n_targets = static_cast<int>(target_cities.size());
    id_paths_no_endpoints = std::vector<std::vector<int>>(n_targets);
    paths_lengths = std::vector<double>(n_targets, -1.0);

    for (int target_city_id = 0; target_city_id < n_targets; ++target_city_id) {
        const auto &target_city = target_cities[target_city_id];

        if (targets_visible[target_city_id]) {
            // The target is directly visible from the source.
            paths_lengths[target_city_id] = source_city.DistanceTo(target_city);
            continue;
        }

        // Add the target to the graph.
        _graph_points[_target_v] = target_city;
        for (auto reflex_vertex_id: targets_visible_reflex_vertices[target_city_id]) {
            AddEdge(reflex_vertex_id, _target_v);
        }

        paths_lengths[target_city_id] = FindShortestPath(id_paths_no_endpoints[target_city_id]);

        // Cleaning after target.
        for (auto reflex_vertex_id: targets_visible_reflex_vertices[target_city_id]) {
            _graph_neigh[reflex_vertex_id].pop_back();
            _graph_costs[reflex_vertex_id].pop_back();
        }
    }
    // Cleaning after source.
    _graph_neigh[_source_v].clear();
    _graph_costs[_source_v].clear();
}

double VisibilityPlanner::FindShortestPath(
    std::vector<int> &id_path_no_endpoints,
    std::optional<int> source_reflex_vertex_id,
    std::optional<int> target_reflex_vertex_id
) {
    int n_graph = _n_graph;
    if (source_reflex_vertex_id) --n_graph;
    if (target_reflex_vertex_id) --n_graph;
    int source_v = source_reflex_vertex_id ? *source_reflex_vertex_id : _source_v;
    int target_v = target_reflex_vertex_id ? *target_reflex_vertex_id : _target_v;

    struct FHNode {
        int id;
        double value;
    };

    struct FHCompare {
        bool operator()(const FHNode &n1, const FHNode &n2) const {
            return n1.value > n2.value;
        }
    };

    using FibonacciHeap = boost::heap::fibonacci_heap<FHNode, boost::heap::compare<FHCompare>>;

    id_path_no_endpoints.clear();

    std::fill(_d.begin(), _d.end(), std::numeric_limits<double>::max());
    std::fill(_c.begin(), _c.end(), Color::kWhite);
    std::fill(_p.begin(), _p.end(), 0);

    _d[source_v] = 0.0;
    _c[source_v] = Color::kGrey;

    FibonacciHeap queue;
    std::vector<FibonacciHeap::handle_type> handles(n_graph);
    handles[source_v] = queue.push({source_v, _graph_points[source_v].DistanceTo(_graph_points[target_v])});

    while (!queue.empty()) {
        const auto &top_node = queue.top();
        int top_v = top_node.id;
        if (top_v == target_v) {
            break;
        }
        queue.pop();
        const auto &costs = _graph_costs[top_v];
        const auto &neigh = _graph_neigh[top_v];
        for (int i = 0; i < neigh.size(); ++i) {
            int v = neigh[i];
            double cost = costs[i];
            if (cost + _d[top_v] >= _d[v]) {
                continue;
            }
            _d[v] = cost + _d[top_v];
            _p[v] = top_v;
            double new_value = _d[v] + _graph_points[v].DistanceTo(_graph_points[target_v]);
            if (_c[v] == Color::kWhite) {
                _c[v] = Color::kGrey;
                handles[v] = queue.push({v, new_value});
                continue;
            }
            if (_c[v] == Color::kBlack) {
                _c[v] = Color::kGrey;
                queue.push({v, new_value});
                continue;
            }
            queue.update(handles[v], {v, new_value});
        }
        _c[top_v] = Color::kBlack;
    }

    // Reconstruct the path.
    id_path_no_endpoints.clear();
    int v = target_v;
    if (target_reflex_vertex_id) {
        id_path_no_endpoints.push_back(v); // add target
    }
    do {
        v = _p[v];
        id_path_no_endpoints.push_back(v);
    } while (v != source_v);
    if (!source_reflex_vertex_id) {
        id_path_no_endpoints.pop_back(); // remove source
    }

    // Reverse the path to get the correct order.
    std::reverse(id_path_no_endpoints.begin(), id_path_no_endpoints.end());

    return _d[target_v];
}

void VisibilityPlanner::Clear() {
    _is_node_convex->clear();
    _mesh_id_to_reflex_id.clear();
    _n_reflex = 0;
    _source_v = 0;
    _target_v = 0;
    _n_graph = 0;
    _graph_neigh.clear();
    _graph_costs.clear();
    _graph_points.clear();
    _d.clear();
    _c.clear();
    _p.clear();
}
