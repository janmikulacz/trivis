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
#include <iomanip>

using namespace trivis_plus;
using namespace trivis_plus::data_loading;

using namespace trivis;
using namespace trivis::geom;

bool data_loading::SaveTriMesh(
    const mesh::TriMesh &mesh,
    const std::string &file,
    std::stringstream *info
) {
    std::ofstream ofs(file);
    if (!ofs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for writing.\n";
        return false;
    }
    ofs << "MESH TRIANGULAR VERSION 1.0\n";
    ofs << "\n[VERTICES]\n";
    ofs << mesh.vertices.size() << "\n";
    for (int ver_id = 0; ver_id < mesh.vertices.size(); ++ver_id) {
        const auto &vertex = mesh.vertices[ver_id];
        ofs << ver_id;
        // point
        ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << vertex.point.x;
        ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << vertex.point.y;
        // edges
        ofs << " " << vertex.edges.size();
        for (int edge_id: vertex.edges) {
            ofs << " " << edge_id;
        }
        // triangles
        ofs << " " << vertex.triangles.size();
        for (int triangle_id: vertex.triangles) {
            ofs << " " << triangle_id;
        }
        // next weakly simple vertex
        ofs << " " << vertex.next_weak_intersect_ver.value_or(-1);
        ofs << "\n";
    }
    ofs << "\n[EDGES]\n";
    ofs << mesh.edges.size() << "\n";
    for (int edge_id = 0; edge_id < mesh.edges.size(); ++edge_id) {
        const auto &edge = mesh.edges[edge_id];
        ofs << edge_id;
        // vertices
        ofs << " " << edge.vertices[0] << " " << edge.vertices[1];
        // triangles
        ofs << " " << edge.triangles.size();
        for (int triangle_id: edge.triangles) {
            ofs << " " << triangle_id;
        }
        // opposites
        ofs << " " << edge.opposites.size();
        for (int ver_id: edge.opposites) {
            ofs << " " << ver_id;
        }
        ofs << "\n";
    }

    ofs << "\n[TRIANGLES]\n";
    ofs << mesh.triangles.size() << "\n";
    for (int triangle_id = 0; triangle_id < mesh.triangles.size(); ++triangle_id) {
        const auto &triangle = mesh.triangles[triangle_id];
        ofs << triangle_id;
        // vertices
        ofs << " " << triangle.vertices[0] << " " << triangle.vertices[1] << " " << triangle.vertices[2];
        // edges
        ofs << " " << triangle.edges[0] << " " << triangle.edges[1] << " " << triangle.edges[2];
        ofs << "\n";
    }
    ofs.close();
    return true;
}

std::optional<mesh::TriMesh> data_loading::LoadTriMesh(
    const std::string &file,
    std::stringstream *info
) {
    std::ifstream ifs(file);
    if (!ifs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for reading.\n";
        return std::nullopt;
    }
    std::string token;
    ifs >> token;
    if (token != "MESH") {
        if (info) *info << "File " << file << " has incorrect format!\n";
        return std::nullopt;
    }
    ifs >> token;
    if (token != "TRIANGULAR") {
        if (info) *info << "File " << file << " has incorrect mesh type " << token << "!\n";
        return std::nullopt;
    }
    ifs >> token;
    if (token != "VERSION") {
        if (info) *info << "File " << file << " has no version info!\n";
        return std::nullopt;
    }
    ifs >> token;
    if (token != "1.0") {
        if (info) *info << "File " << file << " has incorrect version " << token << "!\n";
        return std::nullopt;
    }
    try {
        mesh::TriMesh out_mesh;
        while (!ifs.eof()) {
            ifs >> token;
            if (token == "[NODES]" || token == "[VERTICES]") {
                int n_vertices;
                ifs >> n_vertices;
                out_mesh.vertices.resize(n_vertices);
                for (int i = 0; i < n_vertices; ++i) {
                    int ver_id;
                    ifs >> ver_id;
                    auto &vertex = out_mesh.vertices[ver_id];
                    // point
                    ifs >> vertex.point.x >> vertex.point.y;
                    // edges
                    int size;
                    ifs >> size;
                    vertex.edges.resize(size);
                    for (int &e: vertex.edges) {
                        ifs >> e;
                    }
                    // triangles
                    ifs >> size;
                    vertex.triangles.resize(size);
                    for (int &t: vertex.triangles) {
                        ifs >> t;
                    }
                    // next weakly simple vertex
                    int next_weakly_simple_ver;
                    ifs >> next_weakly_simple_ver;
                    if (next_weakly_simple_ver != -1) {
                        vertex.next_weak_intersect_ver = next_weakly_simple_ver;
                    }
                }

            } else if (token == "[EDGES]") {
                int n_edges;
                ifs >> n_edges;
                out_mesh.edges.resize(n_edges);
                for (int i = 0; i < n_edges; ++i) {
                    int edge_id;
                    ifs >> edge_id;
                    auto &edge = out_mesh.edges[edge_id];
                    // vertices
                    ifs >> edge.vertices[0] >> edge.vertices[1];
                    // triangles
                    int size;
                    ifs >> size;
                    edge.triangles.resize(size);
                    for (int &t: edge.triangles) {
                        ifs >> t;
                    }
                    // opposites
                    ifs >> size;
                    edge.opposites.resize(size);
                    for (int &n: edge.opposites) {
                        ifs >> n;
                    }
                }
            } else if (token == "[TRIANGLES]") {
                int n_triangles;
                ifs >> n_triangles;
                out_mesh.triangles.resize(n_triangles);
                for (int i = 0; i < n_triangles; ++i) {
                    int triangle_id;
                    ifs >> triangle_id;
                    auto &triangle = out_mesh.triangles[triangle_id];
                    // vertices
                    ifs >> triangle.vertices[0] >> triangle.vertices[1] >> triangle.vertices[2];
                    // edges
                    ifs >> triangle.edges[0] >> triangle.edges[1] >> triangle.edges[2];
                }
            } else {
                if (info) *info << "File " << file << ": unknown token " << token << "!\n";
            }
        }
        ifs.close();
        return out_mesh;
    } catch (const std::exception &e) {
        if (info) *info << "Error while reading file " << file << ": " << e.what() << "\n";
        ifs.close();
        return std::nullopt;
    }
}
