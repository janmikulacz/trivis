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


bool SaveTriMesh(
    const trivis::mesh::TriMesh &mesh,
    const std::string &file,
    std::stringstream *info = nullptr
);

std::optional<trivis::mesh::TriMesh> LoadTriMesh(
    const std::string &file,
    std::stringstream *info = nullptr
);

}

#endif //TRIVIS_PLUS_DATA_LOADING_LOAD_MESH_H_
