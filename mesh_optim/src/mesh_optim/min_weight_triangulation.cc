/**
 * File:   min_weight_triangulation.cc
 *
 * Date:   07.07.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "mesh_optim/min_weight_triangulation.h"

#include "trivis/trivis.h"
#include "trivis/utils/simple_clock.h"

#include <list>

using namespace trivis;
using namespace trivis::geom;
using namespace trivis::mesh;
using namespace mesh_optim;

constexpr double kMaxDouble = std::numeric_limits<double>::max();

struct IdVal {
    int id;
    double val;
};

double mesh_optim::ComputeMeshWeight(
    const TriMesh &mesh,
    const std::vector<std::vector<double>> &weights
) {
    double result = 0.0;
    for (const auto &e: mesh.edges) {
        result += weights[e.vertices[0]][e.vertices[1]];
    }
    return result;
}

void AddEdge(
    int id,
    int i,
    int j,
    int k,
    TriTriangle &tri,
    TriMesh &mesh,
    std::vector<std::vector<int>> &edge_ids
) {
    int e_id;
    if (edge_ids[i][j] == -1) {
        // The edge does not exist.
        e_id = static_cast<int>(mesh.edges.size());
        TriEdge e;
        e.vertices = {i, j};
        e.triangles.push_back(static_cast<int>(mesh.triangles.size()));
        e.opposites.push_back(k);
        edge_ids[i][j] = edge_ids[j][i] = e_id;
        mesh.vertices[i].edges.push_back(e_id);
        mesh.vertices[j].edges.push_back(e_id);
        mesh.edges.push_back(std::move(e));
    } else {
        // The edge already exists.
        e_id = edge_ids[i][j];
        auto &e = mesh.edges[e_id];
        e.triangles.push_back(static_cast<int>(mesh.triangles.size()));
        e.opposites.push_back(k);
    }
    tri.edges[id] = e_id;
}

TriMesh RemoveTriangles(
    const TriMesh &mesh,
    const std::vector<bool> &triangles_to_remove,
    std::vector<std::vector<int>> &edge_ids
) {
    TriMesh ret;

    int n_v = static_cast<int>(mesh.vertices.size());
    edge_ids = std::vector<std::vector<int>>(n_v, std::vector<int>(n_v, -1));

    // Copy the vertices.
    ret.vertices.resize(n_v);
    for (int i = 0; i < n_v; ++i) {
        ret.vertices[i].point = mesh.vertices[i].point;
    }

    // Filter triangles.
    int cnt_new = 0; // number of new triangles
    for (int i = 0; i < mesh.triangles.size(); ++i) {
        if (!triangles_to_remove[i]) {
            const auto &tri = mesh.triangles[i];
            TriTriangle tri_new{};
            tri_new.vertices = tri.vertices;
            AddEdge(0, tri_new.vertices[0], tri_new.vertices[1], tri_new.vertices[2], tri_new, ret, edge_ids);
            AddEdge(1, tri_new.vertices[1], tri_new.vertices[2], tri_new.vertices[0], tri_new, ret, edge_ids);
            AddEdge(2, tri_new.vertices[2], tri_new.vertices[0], tri_new.vertices[1], tri_new, ret, edge_ids);
            ret.triangles.push_back(tri_new);
            for (int ver_id: tri_new.vertices) {
                ret.vertices[ver_id].triangles.push_back(cnt_new);
            }
            ++cnt_new;
        }
    }

    return ret;
}

std::vector<int> mesh_optim::GenerateRandomSimplePolygon(
    const TriMesh &mesh,
    int triangle_seed_id,
    std::optional<int> max_size_opt,
    std::vector<bool> &res_poly_used_triangles,
    std::mt19937 &rng
) {

    std::list<int> res_poly;

    struct EdgeInfo {
        int edge, triangle, opposite;
        EdgeInfo(int edge, int triangle, int opposite) : edge(edge), triangle(triangle), opposite(opposite) {};
    };

    const auto &seed_triangle = mesh.triangles[triangle_seed_id];
    res_poly_used_triangles = std::vector<bool>(mesh.triangles.size(), false);
    std::vector<bool> used_vertices(mesh.vertices.size(), false);
    std::vector<std::list<int>::iterator> edge_pos(mesh.edges.size());
    res_poly_used_triangles[triangle_seed_id] = true;
    for (int ver_id: seed_triangle.vertices) {
        used_vertices[ver_id] = true;
        if (mesh.vertices[ver_id].next_weak_intersect_ver) {
            int partner_id = mesh.vertices[ver_id].next_weak_intersect_ver.value();
            while (partner_id != ver_id) {
                used_vertices[partner_id] = true;
                partner_id = mesh.vertices[partner_id].next_weak_intersect_ver.value();
            }
        }
    }

    std::vector<EdgeInfo> edges;
    for (int e_id: seed_triangle.edges) {
        if (mesh.edges[e_id].is_boundary()) {
            continue;
        }
        if (res_poly_used_triangles[mesh.edges[e_id].triangles[0]]) {
            edges.emplace_back(e_id, mesh.edges[e_id].triangles[1], 1);
        } else {
            edges.emplace_back(e_id, mesh.edges[e_id].triangles[0], 0);
        }
    }

    for (int ver_id: seed_triangle.vertices) {
        res_poly.push_back(ver_id);
    }

    while (!edges.empty()) {

        if (max_size_opt && res_poly.size() > max_size_opt) {
            break;
        }

        int rand_e_info_id = std::uniform_int_distribution<int>{0, static_cast<int>(edges.size()) - 1}(rng);
        const auto &rand_e_info = edges[rand_e_info_id];
        int curr_e_id = rand_e_info.edge;
        int curr_tri_id = rand_e_info.triangle;
        int curr_opp_id = mesh.edges[curr_e_id].opposites[rand_e_info.opposite];

        // Remove the edge info from edges.
        edges[rand_e_info_id] = edges[edges.size() - 1];
        edges.pop_back();

        if (res_poly_used_triangles[curr_tri_id] || used_vertices[curr_opp_id]) {
            continue;
        }

        res_poly_used_triangles[curr_tri_id] = true;

        const auto &curr_tri = mesh.triangles[curr_tri_id];
        const auto &curr_e = mesh.edges[curr_e_id];

        if ((res_poly.front() == curr_e.vertices[0] || res_poly.front() == curr_e.vertices[1]) && (res_poly.back() == curr_e.vertices[0] || res_poly.back() == curr_e.vertices[1])) {
            res_poly.emplace_front(curr_opp_id);
        } else {
            for (auto it = res_poly.begin(); it != res_poly.end(); it++) {
                if (*it == curr_e.vertices[0] || *it == curr_e.vertices[1]) {
                    res_poly.emplace(++it, curr_opp_id);
                    break;
                }
            }
        }

        for (int e_id: curr_tri.edges) {
            const auto &e = mesh.edges[e_id];
            if (!e.is_boundary()) {
                if (!res_poly_used_triangles[e.triangles[0]]) {
                    edges.emplace_back(e_id, e.triangles[0], 0);
                }
                if (!res_poly_used_triangles[e.triangles[1]]) {
                    edges.emplace_back(e_id, e.triangles[1], 1);
                }
            }
        }

        for (auto ver_id: curr_tri.vertices) {
            used_vertices[ver_id] = true;
            if (mesh.vertices[ver_id].next_weak_intersect_ver) {
                int partner_id = mesh.vertices[ver_id].next_weak_intersect_ver.value();
                while (partner_id != ver_id) {
                    used_vertices[partner_id] = true;
                    partner_id = mesh.vertices[partner_id].next_weak_intersect_ver.value();
                }
            }
        }
    }

    return {res_poly.begin(), res_poly.end()};
}

// A utility function to find weight of a triangle.
// The weight is considered as perimeter (sum of the lengths of all edges) of the triangle
inline double ComputeWeightTriangle(
    const std::vector<std::vector<double>> &weights,
    int i,
    int j,
    int k
) {
    if (weights[i][j] == kMaxDouble || weights[j][k] == kMaxDouble || weights[k][i] == kMaxDouble) {
        return kMaxDouble;
    }
    return (weights[i][j] + weights[j][k] + weights[k][i]) / 2.0;
}

void AddToTrianglesToMeshRecursive(
    const std::vector<int> &simple_poly_vertices,
    const std::vector<std::vector<IdVal>> &table,
    int i,
    int j,
    std::vector<std::vector<int>> &edge_ids,
    TriMesh &mesh
) {
    if ((i != j) && (table[i][j].val != 0)) {

        int k = table[i][j].id;

        TriTriangle triangle{};
        triangle.vertices[0] = simple_poly_vertices[i];
        triangle.vertices[1] = simple_poly_vertices[k];
        triangle.vertices[2] = simple_poly_vertices[j];

        AddEdge(0, triangle.vertices[0], triangle.vertices[1], triangle.vertices[2], triangle, mesh, edge_ids);
        AddEdge(1, triangle.vertices[1], triangle.vertices[2], triangle.vertices[0], triangle, mesh, edge_ids);
        AddEdge(2, triangle.vertices[2], triangle.vertices[0], triangle.vertices[1], triangle, mesh, edge_ids);

        int n = static_cast<int>(mesh.triangles.size());
        mesh.triangles.push_back(triangle);
        for (int ver_id: triangle.vertices) {
            mesh.vertices[ver_id].triangles.push_back(n);
        }

        AddToTrianglesToMeshRecursive(simple_poly_vertices, table, i, k, edge_ids, mesh);
        AddToTrianglesToMeshRecursive(simple_poly_vertices, table, k, j, edge_ids, mesh);
    }
}

// A Dynamic programming based function to find minimum weight for simple polygon triangulation.
std::vector<std::vector<IdVal>> MinWeightTriangulationOfSimplePolygonDynamicProgramming(
    const std::vector<std::vector<double>> &weights
) {
    int n = static_cast<int>(weights.size());

    // There must be at least 3 points to form a triangle
    if (n < 3) {
        return std::vector<std::vector<IdVal>>{};
    }

    // Table to store results of sub-problems.
    // The value of table[i][j] is the weight of triangulation of points from i to j.
    // The value of table[0][n-1] is the final result.
    std::vector<std::vector<IdVal>> table(n, std::vector<IdVal>(n));

    // Fill table using above recursive formula.
    // Note that the table is filled in diagonal fashion i.e., from diagonal elements to table[0][n-1] which is the result.
    for (int gap = 0; gap < n; ++gap) {
        // std::cout << "gap: " << gap << " / " << n << "\n";
        for (int i = 0, j = gap; j < n; ++i, ++j) {
            if (j < i + 2) {
                table[i][j].val = 0.0;
            } else {
                table[i][j].val = kMaxDouble;
                for (int k = i + 1; k < j; ++k) {
                    double val = table[i][k].val + table[k][j].val + ComputeWeightTriangle(weights, i, j, k);
                    if (val < table[i][j].val) {
                        table[i][j].val = val;
                        table[i][j].id = k;
                    }
                }
            }
        }
    }
    return table;
}

std::vector<std::vector<IdVal>> MinWeightTriangulationOfSimplePolygon(
    const std::vector<int> &simple_poly_vertices,
    const Trivis &simple_poly_vis,
    const std::vector<std::vector<double>> &weights
) {
    int n = static_cast<int>(simple_poly_vertices.size());
    std::vector<int> simple_poly_vertices_inv(simple_poly_vis.mesh().vertices.size(), -1);
    for (int i = 0; i < n; ++i) {
        simple_poly_vertices_inv[simple_poly_vertices[i]] = i;
    }

    std::vector<std::vector<double>> simple_poly_weights(n, std::vector<double>(n, kMaxDouble));
    for (int i = 0; i < simple_poly_vertices.size(); ++i) {
        int ver_id_i = simple_poly_vertices[i];
        std::vector<int> visible_vertices = simple_poly_vis.VisibleVertices(ver_id_i);
        for (int ver_id_j: visible_vertices) {
            int j = simple_poly_vertices_inv[ver_id_j];
            simple_poly_weights[i][j] = simple_poly_weights[j][i] = weights[ver_id_i][ver_id_j];
        }
    }

    auto table = MinWeightTriangulationOfSimplePolygonDynamicProgramming(simple_poly_weights);

    return table;
}

mesh::TriMesh mesh_optim::MinWeightTriangulationIteration(
    const TriMesh &mesh_prev,
    const std::vector<std::vector<double>> &weights,
    std::optional<int> simple_poly_max_size_opt,
    std::mt19937 &rng
) {

    // Generate a random simple polygon based on the previous mesh.
    std::vector<bool> simple_poly_used_triangles;
    int random_tri_id = std::uniform_int_distribution<int>{0, static_cast<int>(mesh_prev.triangles.size()) - 1}(rng);
    std::vector<int> simple_poly_vertices = GenerateRandomSimplePolygon(mesh_prev, random_tri_id, simple_poly_max_size_opt, simple_poly_used_triangles, rng);
    int n = static_cast<int>(simple_poly_vertices.size());

    // Remove the triangles present in the simple polygon from the current mesh.
    std::vector<std::vector<int>> edge_ids, unused;
    auto mesh_new = RemoveTriangles(mesh_prev, simple_poly_used_triangles, edge_ids); // remove used triangles

    // Create a Trivis structure for computing visibility inside the simple polygon.
    simple_poly_used_triangles.flip(); // flip all bool values in the vector
    auto simple_poly_mesh = RemoveTriangles(mesh_prev, simple_poly_used_triangles, unused); // remove unused triangles
    Preorder(simple_poly_mesh, true);
    Trivis simple_poly_vis;
    simple_poly_vis.SetMesh(simple_poly_mesh);

    // Compute the MWT of the simple polygon and add it to the current mesh.
    auto best_triangles_table = MinWeightTriangulationOfSimplePolygon(simple_poly_vertices, simple_poly_vis, weights);
    AddToTrianglesToMeshRecursive(simple_poly_vertices, best_triangles_table, 0, n - 1, edge_ids, mesh_new);
    Preorder(mesh_new);

    return mesh_new;
}

TriMesh mesh_optim::MinWeightTriangulation(
    TriMesh mesh_init,
    const std::vector<std::vector<double>> &weights,
    std::mt19937 &rng,
    std::optional<int> iter_limit_opt,
    std::optional<double> improvement_ratio_limit_opt,
    std::optional<double> time_limit_opt,
    std::optional<int> simple_poly_max_size_opt,
    double *mesh_weight
) {
    utils::SimpleClock clock;

    double orig_weight = ComputeMeshWeight(mesh_init, weights);
    std::cout << "Orig weight: " << orig_weight << std::endl;

    TriMesh mesh_res = std::move(mesh_init);
    // warning: must no longer access mesh_init (it was moved)!

    double prev_weight = orig_weight;
    double new_weight = prev_weight;

    int i = 0;
    while (!iter_limit_opt || iter_limit_opt > i++) {
        if (time_limit_opt && clock.TimeInSeconds() > time_limit_opt) {
            break;
        }
        auto mesh_new = MinWeightTriangulationIteration(mesh_res, weights, simple_poly_max_size_opt, rng);
        new_weight = ComputeMeshWeight(mesh_new, weights);
        double improvement_ratio = (prev_weight - new_weight) / std::abs(prev_weight);
        std::cout << "[" << i << "] New weight: " << new_weight << "; improvement ratio: " << improvement_ratio << ".\n";
        if (improvement_ratio_limit_opt && improvement_ratio < improvement_ratio_limit_opt) {
            break;
        }
        prev_weight = new_weight;
        mesh_res = std::move(mesh_new);
    }
    if (mesh_weight) *mesh_weight = new_weight;

    double percentage_improvement = (orig_weight - new_weight) / std::abs(prev_weight);
    std::cout << "[DONE] New weight: " << new_weight << "; improvement ratio: " << percentage_improvement << ".\n";

    return mesh_res;
}
