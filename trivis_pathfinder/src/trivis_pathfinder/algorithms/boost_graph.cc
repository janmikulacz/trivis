/**
 * File:   boost_graph.cc
 *
 * Date:   21.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_pathfinder/algorithms/boost_graph.h"

using namespace trivis_pathfinder;
using namespace trivis_pathfinder::algorithms;

BoostGraph algorithms::ConstructReflexGraph(
    const std::vector<std::vector<int>> &vis_graph_reflex_reflex,
    const trivis::geom::FPoints &reflex_vertices_points
) {
    int n_reflex_vertices = static_cast<int>(reflex_vertices_points.size());
    BoostGraph graph(n_reflex_vertices);
    std::vector<std::vector<bool>> added(n_reflex_vertices, std::vector<bool>(n_reflex_vertices, false));
    for (int i = 0; i < n_reflex_vertices; ++i) {
        for (int j: vis_graph_reflex_reflex[i]) {
            if (added[i][j] || added[j][i]) {
                continue;
            }
            const auto &pi = reflex_vertices_points[i];
            const auto &pj = reflex_vertices_points[j];
            boost::add_edge(i, j, pi.DistanceTo(pj), graph);
            added[i][j] = true;
            added[j][i] = true;
        }
    }
    return graph;
}
