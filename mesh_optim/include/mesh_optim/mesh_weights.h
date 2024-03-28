/**
 * File:   mesh_weights.h
 *
 * Date:   20.10.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */


#ifndef MESH_OPTIM_MESH_WEIGHTS_H_
#define MESH_OPTIM_MESH_WEIGHTS_H_

#include "trivis/trivis.h"

namespace mesh_optim {

constexpr double kMaxDouble = std::numeric_limits<double>::max();

std::vector<std::vector<double>> ComputeWeightsEdgeLength(
    const trivis::Trivis &vis,
    double eps_dist_diff_collinear
);

std::vector<std::vector<double>> ComputeWeightsEdgeVisibility(
    const trivis::Trivis &vis,
    double max_sample_dist,
    double eps_dist_diff_collinear,
    std::optional<double> vis_radius = std::nullopt,
    bool parallelize = false
);

std::vector<std::vector<double>> ComputeWeightsEdgeVisibilityOptimized(
    const trivis::Trivis &vis,
    double max_sample_dist,
    double eps_dist_diff_collinear,
    double long_edge_penalty,
    std::optional<double> long_edge_threshold_opt = std::nullopt,
    std::optional<int> shortest_edges_node_max_opt = std::nullopt,
    std::optional<double> vis_radius_opt = std::nullopt
);

void NegateWeights(
    std::vector<std::vector<double>> &costs
);

}

#endif //MESH_OPTIM_MESH_WEIGHTS_H_
