/**
 * File:   triangular_meshes.cc
 *
 * Date:   07.07.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "mesh_optim/triangular_meshes.h"

#include <random>

#include "trivis/trivis.h"
#include "trivis/mesh/cdt.h"

#include "mesh_optim/min_weight_triangulation.h"

using namespace trivis;
using namespace trivis::geom;
using namespace mesh_optim;

// A utility function to find weight of a triangle. The weight is considered
// as perimeter (sum of lengths of all edges) of the triangle
inline double ComputeWeightTriangle(
    const std::vector<std::vector<double>> &weights,
    int i,
    int j,
    int k
) {
    return (weights[i][j] + weights[j][k] + weights[k][i]) / 2.0;
}

mesh::TriMesh mesh_optim::GenerateConstrainedDelaunayTriangulation(
    const Trivis &vis
) {
    mesh::TriMesh ret;
    geom::Polygons<double> unused;
    mesh::TriangulateMapCDT(vis.map(), ret, unused);
    return ret;
}

mesh::TriMesh mesh_optim::GenerateMinWeightTriangulation(
    const Trivis &vis,
    const std::vector<std::vector<double>> &weights,
    long random_seed,
    std::optional<int> iter_limit_opt,
    std::optional<double> improvement_ratio_limit_opt,
    std::optional<double> time_limit_opt,
    std::optional<int> simple_poly_max_size_opt,
    double *weight
) {
    auto rng = std::mt19937(random_seed);
    return MinWeightTriangulation(vis.mesh(), weights, rng, iter_limit_opt, improvement_ratio_limit_opt, time_limit_opt, simple_poly_max_size_opt, weight);
}
