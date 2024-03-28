/**
 * File:   min_weight_triangulation.h
 *
 * Date:   07.07.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef MESH_OPTIM_MIN_WEIGHT_TRIANGULATION_H_
#define MESH_OPTIM_MIN_WEIGHT_TRIANGULATION_H_

#include <random>
#include <optional>

#include "trivis/geom/poly_map.h"
#include "trivis/mesh/tri_mesh.h"

namespace mesh_optim {

double ComputeMeshWeight(
    const trivis::mesh::TriMesh &mesh,
    const std::vector<std::vector<double>> &weights
);

std::vector<int> GenerateRandomSimplePolygon(
    const trivis::mesh::TriMesh &mesh,
    int triangle_seed_id,
    std::optional<int> max_size_opt,
    std::vector<bool> &res_poly_used_triangles,
    std::mt19937 &rng
);

trivis::mesh::TriMesh MinWeightTriangulationIteration(
    const trivis::mesh::TriMesh &mesh_prev,
    const std::vector<std::vector<double>> &weights,
    std::optional<int> simple_poly_max_size_opt,
    std::mt19937 &rng
);

trivis::mesh::TriMesh MinWeightTriangulation(
    trivis::mesh::TriMesh mesh_init,
    const std::vector<std::vector<double>> &weights,
    std::mt19937 &rng,
    std::optional<int> iter_limit_opt = 10,
    std::optional<double> improvement_ratio_limit_opt = 1e-4,
    std::optional<double> time_limit_opt = std::nullopt,
    std::optional<int> simple_poly_max_size_opt = std::nullopt,
    double *mesh_weight = nullptr
);

}

#endif //MESH_OPTIM_MIN_WEIGHT_TRIANGULATION_H_
