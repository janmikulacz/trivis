/**
 * File:   triangulations.h
 *
 * Date:   07.07.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef MESH_OPTIM_TRIANGULAR_MESHES_H_
#define MESH_OPTIM_TRIANGULAR_MESHES_H_

#include "trivis/trivis.h"
#include "trivis/mesh/tri_mesh.h"
#include "trivis/geom/poly_map.h"

namespace mesh_optim {

trivis::mesh::TriMesh GenerateConstrainedDelaunayTriangulation(
    const trivis::Trivis &vis
);

trivis::mesh::TriMesh GenerateMinWeightTriangulation(
    const trivis::Trivis &vis,
    const std::vector<std::vector<double>> &weights,
    long random_seed,
    std::optional<int> iter_limit_opt = 10,
    std::optional<double> improvement_ratio_limit_opt = 1e-4,
    std::optional<double> time_limit_opt = std::nullopt,
    std::optional<int> simple_poly_max_size_opt = std::nullopt,
    double *weight = nullptr
);

}

#endif //MESH_OPTIM_TRIANGULAR_MESHES_H_
