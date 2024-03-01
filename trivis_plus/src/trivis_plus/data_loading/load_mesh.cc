/**
 * File:   load_mesh.cc
 *
 * Date:   27.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/data_loading/load_mesh.h"

#include <fstream>
#include <iostream>

bool trivis_plus::data_loading::LoadTriMesh(
    const std::string &file,
    const trivis::geom::FPoints &points,
    trivis::mesh::TriMesh &mesh
) noexcept(false) {
    std::ifstream ifs(file.c_str());
    if (ifs.fail()) {
        return false;
    }
    std::string token;
    while (true) {
        ifs >> token;
        if (ifs.eof()) {
            break;
        }
        if (token == "[NODES]") {
            int n_nodes;
            ifs >> n_nodes;
            mesh.nodes.resize(n_nodes);
            for (int i = 0; i < n_nodes; ++i) {
                int node_id, point_id;
                ifs >> node_id;
                ifs >> point_id;
                auto &node = mesh.nodes[node_id];
                node.point = points[point_id];
                int n_edges;
                ifs >> n_edges;
                node.edges.resize(n_edges);
                for (int j = 0; j < n_edges; ++j) {
                    int edge_id;
                    ifs >> edge_id;
                    node.edges[j] = edge_id;
                }
                int n_triangles;
                ifs >> n_triangles;
                node.triangles.resize(n_triangles);
                for (int j = 0; j < n_triangles; ++j) {
                    int triangle_id;
                    ifs >> triangle_id;
                    node.triangles[j] = triangle_id;
                }
            }
            continue;
        }
        if (token == "[EDGES]") {
            int n_edges;
            ifs >> n_edges;
            mesh.edges.resize(n_edges);
            for (int i = 0; i < n_edges; ++i) {
                int edge_id;
                ifs >> edge_id;
                auto &edge = mesh.edges[edge_id];
                ifs >> edge.nodes[0];
                ifs >> edge.nodes[1];
                int n_triangles;
                ifs >> n_triangles;
                edge.triangles.resize(n_triangles);
                for (int j = 0; j < n_triangles; ++j) {
                    int triangle_id;
                    ifs >> triangle_id;
                    edge.triangles[j] = triangle_id;
                }
                edge.opposites.resize(n_triangles);
                for (int j = 0; j < n_triangles; ++j) {
                    int opposite_id;
                    ifs >> opposite_id;
                    edge.opposites[j] = opposite_id;
                }
            }
            continue;
        }
        if (token == "[TRIANGLES]") {
            int n_triangles;
            ifs >> n_triangles;
            mesh.triangles.resize(n_triangles);
            for (int i = 0; i < n_triangles; ++i) {
                int triangle_id;
                ifs >> triangle_id;
                auto &triangle = mesh.triangles[triangle_id];
                ifs >> triangle.nodes[0];
                ifs >> triangle.nodes[1];
                ifs >> triangle.nodes[2];
                ifs >> triangle.edges[0];
                ifs >> triangle.edges[1];
                ifs >> triangle.edges[2];
            }
            continue;
        }
    }

    return true;
}

std::string trivis_plus::data_loading::LoadTriMeshSafely(
    const std::string &file,
    const trivis::geom::FPoints &points,
    trivis::mesh::TriMesh &mesh
) {
    try {
        if (!LoadTriMesh(file, points, mesh)) {
            return "File " + file + " could not be opened or found.";
        }
    } catch (const std::invalid_argument &e) {
        return "Error while loading " + file + ":\n" + e.what();
    }
    return "ok";
}
