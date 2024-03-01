/**
 * File:   load_mesh.h
 *
 * Date:   29.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DATA_LOADING_LOAD_MESH_H_
#define TRIVIS_PLUS_DATA_LOADING_LOAD_MESH_H_

#include <string>

#include "trivis/geom/geom_types.h"
#include "trivis/mesh/tri_mesh.h"

namespace trivis_plus::data_loading {

bool LoadTriMesh(
    const std::string &file,
    const trivis::geom::FPoints &points,
    trivis::mesh::TriMesh &mesh
) noexcept(false);

std::string LoadTriMeshSafely(
    const std::string &file,
    const trivis::geom::FPoints &points,
    trivis::mesh::TriMesh &mesh
);

}

#endif //TRIVIS_PLUS_DATA_LOADING_LOAD_MESH_H_
