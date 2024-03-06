/**
 * File:   tri_mesh.cc
 *
 * Date:   24.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/mesh/tri_mesh.h"

#include <cassert>
#include <iostream>

#include "trivis/geom/robust_geometry.h"

using namespace trivis;
using namespace trivis::geom;
using namespace trivis::mesh;

void mesh::OrderEdgesInTrianglesCCW(TriMesh &mesh) {
    const auto &edges = mesh.edges;
    auto &triangles = mesh.triangles;
    auto TurnsRight = [mesh](int i, int j, int k) -> bool {
        return geom::TurnsRight(mesh.point(i), mesh.point(j), mesh.point(k));
    };
    // sort edges in each triangle counterclockwise
    for (auto &tri: triangles) {
        int a = tri.edges[0];
        int b = tri.edges[1];
        int a1 = edges[a].vertices[0];
        int a2 = edges[a].vertices[1];
        int b1 = edges[b].vertices[0];
        int b2 = edges[b].vertices[1];
        if ((a1 == b1 && TurnsRight(a2, a1, b2)) ||
            (a1 == b2 && TurnsRight(a2, a1, b1)) ||
            (a2 == b1 && TurnsRight(a1, a2, b2)) ||
            (a2 == b2 && TurnsRight(a1, a2, b1))) {
            tri.edges[0] = b;
            tri.edges[1] = a;
        }
    }
}

void mesh::SplitWeaklyIntersectingNodes(TriMesh &mesh) {
    auto &edges = mesh.edges;
    auto &nodes = mesh.vertices;
    auto &triangles = mesh.triangles;

    int num_nodes = static_cast<int>(nodes.size());
    for (int node_id = 0; node_id < num_nodes; ++node_id) {
        auto &node = nodes[node_id];
        const auto &node_p = node.point;
        int num_obstacle_edges = 0;
        std::vector<int> obstacle_edges_starts;
        for (int edge_id: node.edges) {
            const auto &edge = edges[edge_id];
            if (edge.is_boundary()) {
                const auto &e_other_node_p = nodes[edge.vertices[0] == node_id ? edge.vertices[1] : edge.vertices[0]].point;
                const auto &e_opp_p = nodes[edge.opposites[0]].point;
                if (!TurnsRight(node_p, e_other_node_p, e_opp_p)) {
                    obstacle_edges_starts.push_back(edge_id);
                }
                ++num_obstacle_edges;
            }
        }
        assert(num_obstacle_edges != 0 && num_obstacle_edges % 2 == 0);
        assert(num_obstacle_edges == static_cast<int>(obstacle_edges_starts.size()) * 2);
        int num_new_nodes = num_obstacle_edges / 2;
        if (num_new_nodes < 2) {
            continue;
        }
        std::vector<int> new_nodes_ids;
        for (int i = 0; i < num_new_nodes; ++i) {
            int new_node_id;
            if (i == 0) {
                new_node_id = node_id;
                node.edges.clear();
                node.triangles.clear();
            } else {
                new_node_id = static_cast<int>(nodes.size());
                nodes.emplace_back();
                nodes.back().point = node_p;
            }
            new_nodes_ids.push_back(new_node_id);
            auto &new_node = nodes[new_node_id];
            int curr_e_id = obstacle_edges_starts[i];
            int curr_tri_id = edges[curr_e_id].triangles[0];
            while (true) {
                edges[curr_e_id].vertices[edges[curr_e_id].vertices[0] == node_id ? 0 : 1] = new_node_id; // update edge node
                new_node.edges.push_back(curr_e_id);
                triangles[curr_tri_id].vertices[triangles[curr_tri_id].vertices[0] == node_id ? 0 : (triangles[curr_tri_id].vertices[1] == node_id ? 1 : 2)] = new_node_id;
                new_node.triangles.push_back(curr_tri_id);
                // update current edge
                const auto &curr_tri = triangles[curr_tri_id];
                int curr_tri_curr_e_idx = curr_tri.edges[0] == curr_e_id ? 0 : (curr_tri.edges[1] == curr_e_id ? 1 : 2);
                int next_e_id = curr_tri.edges[(curr_tri_curr_e_idx + 1) % 3];
                edges[next_e_id].opposites[edges[next_e_id].opposites[0] == node_id ? 0 : 1] = new_node_id;
                curr_e_id = curr_tri.edges[(curr_tri_curr_e_idx + 2) % 3];
                if (!edges[curr_e_id].is_boundary()) {
                    // update current triangle
                    curr_tri_id = edges[curr_e_id].triangles[0] == curr_tri_id ? edges[curr_e_id].triangles[1] : edges[curr_e_id].triangles[0];
                } else {
                    // reached the last edge
                    edges[curr_e_id].vertices[edges[curr_e_id].vertices[0] == node_id ? 0 : 1] = new_node_id;
                    new_node.edges.push_back(curr_e_id);
                    break;
                }
            }
        }
        // Save split partners.
        for (int i_prev = static_cast<int>(new_nodes_ids.size()) - 1, i = 0; i < new_nodes_ids.size(); i_prev = i++) {
            nodes[new_nodes_ids[i_prev]].next_wi_vertex = new_nodes_ids[i];
        }
    }
}

void mesh::OrderEdgesAndTrianglesInNodesCCW(TriMesh &mesh) {
    const auto &edges = mesh.edges;
    const auto &triangles = mesh.triangles;
    auto &nodes = mesh.vertices;

    // sort edges and triangles in each node counterclockwise
    for (int node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        auto &node = nodes[node_idx];
        const auto &node_p = node.point;
        assert(node.edges.empty() == node.triangles.empty());
        if (node.edges.empty()) {
            continue;
        }
        std::cerr.precision(std::numeric_limits<double>::max_digits10);
        std::vector<int> new_edges;
        std::vector<int> new_triangles;
        int curr_e_idx = -1;
        for (int e_idx: node.edges) {
            const auto &edge = edges[e_idx];
            assert(edge.vertices[0] == node_idx || edge.vertices[1] == node_idx);
            if (edge.is_boundary()) {
                const auto &e_other_node_p = nodes[edge.vertices[0] == node_idx ? edge.vertices[1] : edge.vertices[0]].point;
                const auto &e_opp_p = nodes[edge.opposites[0]].point;
                if (!TurnsRight(node_p, e_other_node_p, e_opp_p)) {
                    curr_e_idx = e_idx;
                    break;
                }
            }
        }
        assert(curr_e_idx != -1);
        int curr_tri_idx = edges[curr_e_idx].triangles[0];
        while (true) {
            new_edges.push_back(curr_e_idx);
            new_triangles.push_back(curr_tri_idx);
            // update current edge
            const auto &curr_tri = triangles[curr_tri_idx];
            int curr_tri_curr_e_idx = curr_tri.edges[0] == curr_e_idx ? 0 : (curr_tri.edges[1] == curr_e_idx ? 1 : 2);
            curr_e_idx = curr_tri.edges[(curr_tri_curr_e_idx + 2) % 3];
            if (!edges[curr_e_idx].is_boundary()) {
                // update current triangle
                curr_tri_idx = edges[curr_e_idx].triangles[0] == curr_tri_idx ? edges[curr_e_idx].triangles[1] : edges[curr_e_idx].triangles[0];
            } else {
                // reached the last edge
                new_edges.push_back(curr_e_idx);
                break;
            }
        }
        assert(new_edges.size() == node.edges.size());
        assert(new_triangles.size() == node.triangles.size());
        node.edges = std::move(new_edges);
        node.triangles = std::move(new_triangles);
    }
}

geom::FPolygons mesh::Mesh2Polygons(const TriMesh &mesh) {
    geom::FPolygons ret;
    for (const auto &triangle: mesh.triangles) {
        geom::FPolygon polygon;
        for (int node_id: triangle.vertices) {
            polygon.push_back(mesh.point(node_id));
        }
        ret.push_back(polygon);
    }
    return ret;
}