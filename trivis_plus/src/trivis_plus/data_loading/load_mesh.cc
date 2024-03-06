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

bool trivis_plus::data_loading::SaveTriMesh(
    const trivis::mesh::TriMesh &mesh,
    const std::string &file,
    std::stringstream *info
) {
    std::fstream fs;
    fs.open(file, std::ios::out | std::ios::trunc);
    if (!fs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for writing.\n";
        return false;
    }
    fs << "MESH TRIANGULAR VERSION 1.0\n";
    fs << "\n[vertexS]\n";
    fs << mesh.vertices.size() << "\n";
    for (int ver_id = 0; ver_id < mesh.vertices.size(); ++ver_id) {
        const auto &vertex = mesh.vertices[ver_id];
        fs << ver_id;
        // point
        fs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << vertex.point.x;
        fs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << vertex.point.y;
        // edges
        fs << " " << vertex.edges.size();
        for (int edge_id: vertex.edges) {
            fs << " " << edge_id;
        }
        // triangles
        fs << " " << vertex.triangles.size();
        for (int triangle_id: vertex.triangles) {
            fs << " " << triangle_id;
        }
        // split partner
        fs << " " << vertex.next_wi_vertex.value_or(-1);
        fs << "\n";
    }
    fs << "\n[EDGES]\n";
    fs << mesh.edges.size() << "\n";
    for (int edge_id = 0; edge_id < mesh.edges.size(); ++edge_id) {
        const auto &edge = mesh.edges[edge_id];
        fs << edge_id;
        // vertices
        fs << " " << edge.vertices[0] << " " << edge.vertices[1];
        // triangles
        fs << " " << edge.triangles.size();
        for (int triangle_id: edge.triangles) {
            fs << " " << triangle_id;
        }
        // opposites
        fs << " " << edge.opposites.size();
        for (int ver_id: edge.opposites) {
            fs << " " << ver_id;
        }
        fs << "\n";
    }

    fs << "\n[TRIANGLES]\n";
    fs << mesh.triangles.size() << "\n";
    for (int triangle_id = 0; triangle_id < mesh.triangles.size(); ++triangle_id) {
        const auto &triangle = mesh.triangles[triangle_id];
        fs << triangle_id;
        // vertices
        fs << " " << triangle.vertices[0] << " " << triangle.vertices[1] << " " << triangle.vertices[2];
        // edges
        fs << " " << triangle.edges[0] << " " << triangle.edges[1] << " " << triangle.edges[2];
        fs << "\n";
    }
    fs.close();
    return true;
}

std::optional<trivis::mesh::TriMesh> trivis_plus::data_loading::LoadTriMesh(
    const std::string &file,
    std::stringstream *info
) {
    std::fstream fs;
    fs.open(file, std::ios::in);
    if (!fs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for reading.\n";
        return std::nullopt;
    }
    std::string token;
    fs >> token;
    if (token != "MESH") {
        if (info) *info << "File " << file << " has incorrect format!\n";
        return std::nullopt;
    }
    fs >> token;
    if (token != "TRIANGULAR") {
        if (info) *info << "File " << file << " has incorrect mesh type " << token << "!\n";
        return std::nullopt;
    }
    fs >> token;
    if (token != "VERSION") {
        if (info) *info << "File " << file << " has no version info!\n";
        return std::nullopt;
    }
    fs >> token;
    if (token != "1.0") {
        if (info) *info << "File " << file << " has incorrect version " << token << "!\n";
        return std::nullopt;
    }
    try {
        trivis::mesh::TriMesh out_mesh;
        while (!fs.eof()) {
            fs >> token;
            if (token == "[vertexS]") {
                int n_vertices;
                fs >> n_vertices;
                out_mesh.vertices.resize(n_vertices);
                for (int i = 0; i < n_vertices; ++i) {
                    int ver_id;
                    fs >> ver_id;
                    auto &vertex = out_mesh.vertices[ver_id];
                    // point
                    fs >> vertex.point.x >> vertex.point.y;
                    // edges
                    int size;
                    fs >> size;
                    vertex.edges.resize(size);
                    for (int &e: vertex.edges) {
                        fs >> e;
                    }
                    // triangles
                    fs >> size;
                    vertex.triangles.resize(size);
                    for (int &t: vertex.triangles) {
                        fs >> t;
                    }
                    // split partner
                    int next_weakly_simple_node;
                    fs >> next_weakly_simple_node;
                    if (next_weakly_simple_node != -1) {
                        vertex.next_wi_vertex = next_weakly_simple_node;
                    }
                }

            } else if (token == "[EDGES]") {
                int n_edges;
                fs >> n_edges;
                out_mesh.edges.resize(n_edges);
                for (int i = 0; i < n_edges; ++i) {
                    int edge_id;
                    fs >> edge_id;
                    auto &edge = out_mesh.edges[edge_id];
                    // vertices
                    fs >> edge.vertices[0] >> edge.vertices[1];
                    // triangles
                    int size;
                    fs >> size;
                    edge.triangles.resize(size);
                    for (int &t: edge.triangles) {
                        fs >> t;
                    }
                    // opposites
                    fs >> size;
                    edge.opposites.resize(size);
                    for (int &n: edge.opposites) {
                        fs >> n;
                    }
                }
            } else if (token == "[TRIANGLES]") {
                int n_triangles;
                fs >> n_triangles;
                out_mesh.triangles.resize(n_triangles);
                for (int i = 0; i < n_triangles; ++i) {
                    int triangle_id;
                    fs >> triangle_id;
                    auto &triangle = out_mesh.triangles[triangle_id];
                    // vertices
                    fs >> triangle.vertices[0] >> triangle.vertices[1] >> triangle.vertices[2];
                    // edges
                    fs >> triangle.edges[0] >> triangle.edges[1] >> triangle.edges[2];
                }
            } else {
                if (info) *info << "File " << file << ": unknown token " << token << "!\n";
            }
        }
        fs.close();
        return out_mesh;
    } catch (const std::exception &e) {
        if (info) *info << "Error while reading file " << file << ": " << e.what() << "\n";
        fs.close();
        return std::nullopt;
    }
}
