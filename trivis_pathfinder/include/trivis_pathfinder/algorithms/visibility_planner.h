/**
 * File:   visibility_planner.h
 *
 * Date:   24.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PATHFINDER_VISIBILITY_PLANNER_H_
#define TRIVIS_PATHFINDER_VISIBILITY_PLANNER_H_

#include "trivis_pathfinder/utils/status.h"
#include "trivis/trivis.h"

namespace trivis_pathfinder::algorithms {

class VisibilityPlanner {
public:

    void Init(
        const std::vector<bool> &is_node_convex,
        const std::vector<int> &mesh_id_to_reflex_id,
        const trivis::geom::FPoints &reflex_points,
        const std::vector<std::vector<int>> &vis_graph_reflex_reflex
    );

    double FindShortestPath(
        int source_reflex_vertex_id,
        int target_reflex_vertex_id,
        std::vector<int> &id_path
    );

    utils::StatusWithResult<double> FindShortestPath(
        const trivis::geom::FPoint &source_city,
        const trivis::geom::FPoint &target_city,
        const trivis::Trivis &vis,
        std::vector<int> &id_path_no_endpoints
    );

    utils::StatusWithResult<double> FindShortestPath(
        const trivis::geom::FPoint &source_city,
        const trivis::geom::FPoint &target_city,
        bool target_visible,
        const std::vector<int> &source_visible_reflex_vertices,
        const std::vector<int> &target_visible_reflex_vertices,
        std::vector<int> &id_path_no_endpoints
    );

    utils::Status FindShortestPaths(
        const trivis::geom::FPoint &source_city,
        const trivis::geom::FPoints &target_cities,
        const trivis::Trivis &vis,
        std::vector<std::vector<int>> &id_paths_no_endpoints,
        std::vector<double> &paths_lengths
    );

    utils::Status FindShortestPaths(
        const trivis::geom::FPoint &source_city,
        const trivis::geom::FPoints &target_cities,
        std::vector<bool> &targets_visible,
        const std::vector<int> &source_visible_reflex_vertices,
        const std::vector<std::vector<int>> &targets_visible_reflex_vertices,
        std::vector<std::vector<int>> &id_paths_no_endpoints,
        std::vector<double> &paths_lengths
    );

    void Clear();

private:

    enum class Color : int {
        kWhite = 0,
        kGrey = 1,
        kBlack = 2,
    };

    std::optional<std::vector<bool>> _is_node_convex;
    std::vector<int> _mesh_id_to_reflex_id;
    int _n_reflex = 0;
    int _source_v = 0;
    int _target_v = 0;
    int _n_graph = 0;
    std::vector<std::vector<int>> _graph_neigh;
    std::vector<std::vector<double>> _graph_costs;
    trivis::geom::FPoints _graph_points;
    std::vector<double> _d;
    std::vector<Color> _c;
    std::vector<int> _p;

    void AddEdge(
        int source_graph_v,
        int target_graph_v
    ) {
        _graph_neigh[source_graph_v].push_back(target_graph_v);
        _graph_costs[source_graph_v].push_back(_graph_points[source_graph_v].DistanceTo(_graph_points[target_graph_v]));
    }

    double FindShortestPath(
        std::vector<int> &id_path_no_endpoints,
        std::optional<int> source_reflex_vertex_id = std::nullopt,
        std::optional<int> target_reflex_vertex_id = std::nullopt
    );
};

}

#endif //TRIVIS_PATHFINDER_VISIBILITY_PLANNER_H_
