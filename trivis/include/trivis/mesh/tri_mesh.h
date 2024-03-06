/**
 * File:   tri_mesh.h
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_MESH_TRI_MESH_H_
#define TRIVIS_MESH_TRI_MESH_H_

#include <array>
#include <optional>

#include "trivis/geom/geom_types.h"

namespace trivis::mesh {

struct TriVertex {
    std::optional<int> next_wi_vertex = std::nullopt;
    geom::FPoint point;
    std::vector<int> edges;
    std::vector<int> triangles;
};

struct TriEdge {
    std::array<int, 2> vertices;
    std::vector<int> triangles;
    std::vector<int> opposites;

    [[nodiscard]] bool is_boundary() const { return triangles.size() == 1; }
};

struct TriTriangle {
    std::array<int, 3> edges;
    std::array<int, 3> vertices;
};

struct TriMesh {
    std::vector<TriVertex> vertices;
    std::vector<TriEdge> edges;
    std::vector<TriTriangle> triangles;

    [[nodiscard]] const auto &point(int id) const { return vertices[id].point; }
};

[[nodiscard]] inline int OppositeVertex(
    const TriMesh &mesh,
    int tri_id,
    int edge_id
) {
    return mesh.triangles[tri_id].vertices[0]
           + mesh.triangles[tri_id].vertices[1]
           + mesh.triangles[tri_id].vertices[2]
           - mesh.edges[edge_id].vertices[0]
           - mesh.edges[edge_id].vertices[1];
}

void OrderEdgesInTrianglesCCW(TriMesh &mesh);

/**
 * -------------
 * |\#########/|
 * | \###O1##/ |
 * |  \#####/  |
 * |   \###/   |
 * |    \#/    |
 * | T1  X  T2 |
 * |    /#\    |
 * |   /###\   |
 * |  /#####\  |
 * | /###O2##\ |
 * |/#########\|
 * -------------
 * O1, O2 ... obstacles
 * T1, T2 ... triangles, part of the triangulation
 * X ... weak self-intersection node shared by both triangles
 *
 * In theory, looking from T1 to T2 through X should create a 1-dimensional antenna.
 * However, this is usually undesired, and TriVis is not designed to handle these cases properly.
 * Therefore, this function finds all weak self-intersection vertices in the mesh and splits them into X1 and X2.
 * X1 belongs to T1, and X2 belongs to T2.
 * Although X1 and X2 have the same coordinates, they are not visible to each other.
 * This function works even when the weak self-intersection node is shared between more than two triangles.
 *
 */
void SplitWeaklyIntersectingNodes(TriMesh &mesh);

void OrderEdgesAndTrianglesInNodesCCW(TriMesh & mesh);

inline void Preorder(
    TriMesh &mesh,
    bool is_strictly_simple = false
) {
    OrderEdgesInTrianglesCCW(mesh);
    if (!is_strictly_simple) {
        SplitWeaklyIntersectingNodes(mesh);
    }
    OrderEdgesAndTrianglesInNodesCCW(mesh);
}

geom::FPolygons Mesh2Polygons(const TriMesh &mesh);

inline bool operator==(const TriVertex &v0, const TriVertex &v1) {
    return v0.point == v1.point && v0.edges == v1.edges && v0.triangles == v1.triangles && v0.next_wi_vertex == v1.next_wi_vertex;
}

inline bool operator!=(const TriVertex &v0, const TriVertex &v1) {
    return !(v0 == v1);
}

inline bool operator==(const TriEdge &e0, const TriEdge &e1) {
    return e0.vertices == e1.vertices && e0.triangles == e1.triangles && e0.opposites == e1.opposites;
}

inline bool operator!=(const TriEdge &e0, const TriEdge &e1) {
    return !(e0 == e1);
}

inline bool operator==(const TriTriangle &t0, const TriTriangle &t1) {
    return t0.vertices == t1.vertices && t0.edges == t1.edges;
}

inline bool operator!=(const TriTriangle &t0, const TriTriangle &t1) {
    return !(t0 == t1);
}

inline bool operator==(const TriMesh &m0, const TriMesh &m1) {
    return m0.vertices == m1.vertices && m0.edges == m1.edges && m0.triangles == m1.triangles;
}

inline bool operator!=(const TriMesh &m0, const TriMesh &m1) {
    return !(m0 == m1);
}

}

#endif //TRIVIS_MESH_TRI_MESH_H_
