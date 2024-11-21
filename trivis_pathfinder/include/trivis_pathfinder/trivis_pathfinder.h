/**
 * File:   trivis_pathfinder.h
 *
 * Date:   27.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PATHFINDER_PATH_FINDER_H_
#define TRIVIS_PATHFINDER_PATH_FINDER_H_

#include "trivis_pathfinder/status.h"
#include "trivis_pathfinder/visibility_planner.h"
#include "trivis_pathfinder/boost_graph.h"

#include "trivis/trivis.h"

namespace trivis_pathfinder {

class TrivisPathfinder {

public:

    Status ConstructReflexVisibilityGraph(const trivis::Trivis &vis);

    Status PrecomputeReflexShortestPaths();

    Status ConstructCitiesVisibilityGraph(
        const trivis::Trivis &vis,
        trivis::geom::FPoints cities
    );

    Status PrecomputeCitiesShortestPaths();

    StatusWithResult<double> ShortestPathReflex(
        const trivis::Trivis &vis,
        int source_reflex_id,
        int target_reflex_id,
        trivis::geom::FPoints *point_path = nullptr,
        std::vector<int> *id_path_no_endpoints = nullptr
    );

    StatusWithResult<double> ShortestPathCities(
        const trivis::Trivis &vis,
        int source_city_id,
        int target_city_id,
        trivis::geom::FPoints *point_path = nullptr,
        std::vector<int> *id_path_no_endpoints = nullptr
    );

    StatusWithResult<double> ShortestPathPoints(
        const trivis::Trivis &vis,
        const trivis::geom::FPoint &source_point,
        const trivis::geom::FPoint &target_point,
        trivis::geom::FPoints *point_path = nullptr,
        std::vector<int> *id_path_no_endpoints = nullptr
    );

    [[nodiscard]] const auto &VisibleReflexFromReflex(int reflex_id) {
        return _vis_graph_reflex_reflex[reflex_id];
    }

    [[nodiscard]] const auto &VisibleCitiesFromReflex(int reflex_id) {
        return _vis_graph_reflex_city[reflex_id];
    }

    [[nodiscard]] const auto &VisibleReflexFromCity(int city_id) {
        return _vis_graph_city_reflex[city_id];
    }

    [[nodiscard]] const auto &VisibleCitiesFromCity(int city_id) {
        return _vis_graph_city_city[city_id];
    }

    [[nodiscard]] bool IsVisibleReflexReflex(int reflex_id_i, int reflex_id_j) {
        return _vis_graph_reflex_reflex[reflex_id_i][reflex_id_j];
    }

    [[nodiscard]] bool IsVisibleReflexCity(int reflex_id, int city_id) {
        return _vis_graph_reflex_city[reflex_id][city_id];
    }

    [[nodiscard]] bool IsVisibleCityReflex(int city_id, int reflex_id) {
        return _vis_graph_city_reflex[city_id][reflex_id];
    }

    [[nodiscard]] bool IsVisibleCityCity(int city_id_i, int city_id_j) {
        return _vis_graph_city_city[city_id_i][city_id_j];
    }

    [[nodiscard]] bool has_vis_graph_reflex() const {
        return _has_vis_graph_reflex;
    }

    [[nodiscard]] int n_reflex() const {
        return _n_reflex;
    }

    [[nodiscard]] const auto &is_node_convex() const {
        return _is_node_convex;
    }

    [[nodiscard]] const auto &reflex_points() const {
        return _reflex_points;
    }

    [[nodiscard]] const auto &reflex_id_to_mesh_id() const {
        return _reflex_id_to_mesh_id;
    }

    [[nodiscard]] const auto &mesh_id_to_reflex_id() const {
        return _mesh_id_to_reflex_id;
    }

    [[nodiscard]] bool has_shortest_paths_reflex() const {
        return _has_shortest_paths_reflex;
    }

    [[nodiscard]] bool has_vis_graph_cities() const {
        return _has_vis_graph_cities;
    }

    [[nodiscard]] int n_cities() const {
        return _n_cities;
    }

    [[nodiscard]] const auto &cities() const {
        return _cities;
    }

    [[nodiscard]] bool has_shortest_paths_cities() const {
        return _has_shortest_paths_cities;
    }

    void Clear();

private:

    /// ##############################################
    /// ## FILLED BY ConstructReflexVisibilityGraph ##
    /// ##############################################

    // True if the below variables were initialized/filled, i.e., ConstructReflexVisibilityGraph was called.
    bool _has_vis_graph_reflex = false;

    /// === REFLEX VERTICES ===

    // The number of reflex vertices. The reflex vertices are indexed from 0 to _n_reflex - 1.
    int _n_reflex = 0;
    // _is_node_convex[i] is true iff the i-th mesh node is convex.
    std::optional<std::vector<bool>> _is_node_convex;
    // A vector of points that correspond to reflex vertices.
    trivis::geom::FPoints _reflex_points;
    // A mapping from reflex vertex ids to mesh node ids.
    std::vector<int> _reflex_id_to_mesh_id;
    // A mapping from mesh node ids to reflex vertex ids. If the i-th mesh node is convex, then _mesh_id_to_reflex_id[i] == -1.
    std::vector<int> _mesh_id_to_reflex_id;

    /// === REFLEX-REFLEX VISIBILITY GRAPH ===

    // _vis_graph_reflex_reflex[i] contains j iff _reflex_points[j] is visible from _reflex_points[i].
    std::vector<std::vector<int>> _vis_graph_reflex_reflex;
    // _vis_graph_reflex_reflex_bool_matrix[i][j] == true iff _reflex_points[j] is visible from _reflex_points[i].
    std::vector<std::vector<bool>> _vis_graph_reflex_reflex_bool_matrix;

    /// === PLANNING CLASSES ===

    // Structure used by boost's Dijkstra algorithm.
    trivis_pathfinder::BoostGraph _reflex_graph;
    // The visibility planner.
    trivis_pathfinder::VisibilityPlanner _vis_planner;

    /// #############################################
    /// ## FILLED BY PrecomputeReflexShortestPaths ##
    /// #############################################

    // True if the below variables were initialized/filled, i.e., PrecomputeReflexShortestPaths was called.
    bool _has_shortest_paths_reflex = false;

    /// === REFLEX-REFLEX SHORTEST PATHS ===

    // _shortest_paths_lengths_reflex[i][j] is the length of the shortest path between _reflex_points[i] and _reflex_points[j].
    std::vector<std::vector<double>> _shortest_paths_lengths_reflex;
    // _shortest_paths_predecessors_reflex[i] is the predecessor structure of the shortest paths from _reflex_points[i] to all other reflex points.
    std::vector<std::vector<int>> _shortest_paths_predecessors_reflex;

    /// ##############################################
    /// ## FILLED BY ConstructReflexVisibilityGraph ##
    /// ##############################################

    // True if the below variables were initialized/filled, i.e., ConstructReflexVisibilityGraph was called.
    bool _has_vis_graph_cities = false;

    // The number of cities. Cities are indexed from 0 to _n_cities - 1.
    int _n_cities = 0;
    // A vector of points that correspond to reflex vertices.
    trivis::geom::FPoints _cities;
    // A mapping from city ids to the mesh.
    std::vector<std::optional<trivis::Trivis::PointLocationResult>> _cities_pl_results;

    /// === CITY-REFLEX / REFLEX-CITY VISIBILITY GRAPH ===

    // _vis_graph_city_reflex[i] contains j iff _reflex_points[j] is visible from _cites[i].
    std::vector<std::vector<int>> _vis_graph_city_reflex;
    // _vis_graph_city_reflex_bool_matrix[i][j] == true iff _reflex_points[j] is visible from _cities[i].
    std::vector<std::vector<bool>> _vis_graph_city_reflex_bool_matrix;
    // _vis_graph_reflex_city[i] contains j iff _cities[j] is visible from _reflex_points[i].
    std::vector<std::vector<int>> _vis_graph_reflex_city;
    // _vis_graph_reflex_city_bool_matrix[i][j] == true iff _cities[j] is visible from _reflex_points[i].
    std::vector<std::vector<bool>> _vis_graph_reflex_city_bool_matrix;

    /// === CITY-CITY VISIBILITY GRAPH ===

    // _vis_graph_city_city[i] contains j iff _cites[j] is visible from _cites[i].
    std::vector<std::vector<int>> _vis_graph_city_city;
    // _vis_graph_city_city_bool_matrix[i][j] == true iff _cities[j] is visible from _cities[i].
    std::vector<std::vector<bool>> _vis_graph_city_city_bool_matrix;

    /// #############################################
    /// ## FILLED BY PrecomputeCitiesShortestPaths ##
    /// #############################################

    // True if the below variables were initialized/filled, i.e., PrecomputeCitiesShortestPaths was called.
    bool _has_shortest_paths_cities = false;

    /// === CITY-CITY SHORTEST PATHS ===

    // _shortest_paths_lengths_cities[i][j] is the length of the shortest path between _cities[i] and _cities[j].
    std::vector<std::vector<double>> _shortest_paths_lengths_cities;
    // _shortest_paths_predecessors_cities[i] is the predecessor structure of the shortest paths from _cities[i] to all other cities.
    std::vector<std::vector<int>> _shortest_paths_predecessors_cities;

};

}

#endif //TRIVIS_PATHFINDER_PATH_FINDER_H_
