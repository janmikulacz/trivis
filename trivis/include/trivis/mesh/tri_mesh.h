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

#include "trivis/geom/geom_types.h"

namespace trivis::mesh {

struct TriNode {
    int next_split_partner = -1; // invalid
    geom::FPoint point;
    std::vector<int> edges;
    std::vector<int> triangles;
};

struct TriEdge {
    std::array<int, 2> nodes;
    std::vector<int> triangles;
    std::vector<int> opposites;
    [[nodiscard]] bool is_obstacle() const { return triangles.size() == 1; }
};

struct TriTriangle {
    std::array<int, 3> edges;
    std::array<int, 3> nodes;
};

struct TriMesh {
    std::vector<TriNode> nodes;
    std::vector<TriEdge> edges;
    std::vector<TriTriangle> triangles;
    [[nodiscard]] const auto &point(int id) const { return nodes[id].point; }
};

[[nodiscard]] inline int OppositeNode(
    const TriMesh &mesh,
    int tri_id,
    int e_id
) {
    return mesh.triangles[tri_id].nodes[0] + mesh.triangles[tri_id].nodes[1] + mesh.triangles[tri_id].nodes[2] - mesh.edges[e_id].nodes[0] - mesh.edges[e_id].nodes[1];
}

void SortEdgesInTriangles(TriMesh &mesh);

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
 * Therefore, this function finds all weak self-intersection nodes in the mesh and splits them into X1 and X2.
 * X1 belongs to T1, and X2 belongs to T2.
 * Although X1 and X2 have the same coordinates, they are not visible to each other.
 * This function works even when the weak self-intersection node is shared between more than two triangles.
 *
 */
void SplitWeakSelfIntersectionNodes(TriMesh &mesh);

void SortEdgesAndTrianglesInNodes(TriMesh &mesh);

inline void Reorganize(TriMesh &mesh, bool is_strictly_simple = false) {
    SortEdgesInTriangles(mesh);
    if (!is_strictly_simple) {
        SplitWeakSelfIntersectionNodes(mesh);
    }
    SortEdgesAndTrianglesInNodes(mesh);
}

geom::FPolygons Mesh2Polygons(const TriMesh &mesh);

inline bool operator==(const TriNode &n0, const TriNode &n1) {
    return n0.point == n1.point && n0.edges == n1.edges && n0.triangles == n1.triangles && n0.next_split_partner == n1.next_split_partner;
}

inline bool operator!=(const TriNode &n0, const TriNode &n1) {
    return !(n0 == n1);
}

inline bool operator==(const TriEdge &e0, const TriEdge &e1) {
    return e0.nodes == e1.nodes && e0.triangles == e1.triangles && e0.opposites == e1.opposites;
}

inline bool operator!=(const TriEdge &e0, const TriEdge &e1) {
    return !(e0 == e1);
}

inline bool operator==(const TriTriangle &t0, const TriTriangle &t1) {
    return t0.nodes == t1.nodes && t0.edges == t1.edges;
}

inline bool operator!=(const TriTriangle &t0, const TriTriangle &t1) {
    return !(t0 == t1);
}

inline bool operator==(const TriMesh &m0, const TriMesh &m1) {
    return m0.nodes == m1.nodes && m0.edges == m1.edges && m0.triangles == m1.triangles;
}

inline bool operator!=(const TriMesh &m0, const TriMesh &m1) {
    return !(m0 == m1);
}

}

#endif //TRIVIS_MESH_TRI_MESH_H_
