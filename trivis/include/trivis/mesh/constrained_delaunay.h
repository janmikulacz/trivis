/**
 * File:   constrained_delaunay.h
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_MESH_CONSTRAINED_DELAUNAY_H_
#define TRIVIS_MESH_CONSTRAINED_DELAUNAY_H_

#include "trivis/geom/poly_map.h"

#include "trivis/mesh/tri_mesh.h"

namespace trivis::mesh {

void TriangulateMapConstrainedDelaunay(
    const geom::PolyMap &map,
    TriMesh &mesh,
    geom::FPolygons &triangles
);

void TriangulateMapConstrainedDelaunay(
    const geom::FPolygons &polygons,
    geom::FPolygons &triangles
);

void TriangulateMapConstrainedConformingDelaunay(
    const geom::PolyMap &map,
    double max_circumradius,
    geom::FPolygons &triangles
);

}

#endif //TRIVIS_MESH_CONSTRAINED_DELAUNAY_H_
