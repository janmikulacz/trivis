/**
 * File:   cdt.h
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_MESH_CDT_H_
#define TRIVIS_MESH_CDT_H_

#include "trivis/geom/poly_map.h"

#include "trivis/mesh/tri_mesh.h"

namespace trivis::mesh {

void TriangulateMapCDT(
    const geom::PolyMap &map,
    TriMesh &mesh,
    geom::FPolygons &triangles
);

void TriangulateMapCDT(
    const geom::FPolygons &polygons,
    geom::FPolygons &triangles
);

void TriangulateMapCCDT(
    const geom::PolyMap &map,
    double max_circumradius,
    geom::FPolygons &triangles
);

}

#endif //TRIVIS_MESH_CDT_H_
