/**
 * File:   trivis.cc
 *
 * Date:   29.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/trivis.h"

#include <fstream>
#include <iomanip>
#include <cassert>

#include "trivis/mesh/cdt.h"

#include "trivis/geom/generic_geom_utils.h"
#include "trivis/geom/robust_geometry.h"
#include "trivis/geom/intersections.h"

using namespace trivis;
using namespace trivis::geom;

///================///
/// PUBLIC MEMBERS ///
///================///

/// ##### INITIALIZATION ##### ///

void Trivis::SetMap(PolyMap map) {
    _mesh.nodes.clear();
    _mesh.edges.clear();
    _mesh.triangles.clear();
    _triangles.clear();
    _has_mesh = false;
    _limits = map.limits();
    _map = std::move(map);
    _has_map = true;
    _has_pl = false;
}

void Trivis::ConstructMeshCDT() {
    assert(_has_map);
    mesh::TriangulateMapCDT(_map, _mesh, _triangles);
    _has_mesh = true;
}

void Trivis::SetMesh(
    mesh::TriMesh mesh,
    std::optional<geom::PolyMap> map
) {
    _triangles.clear();
    _triangles.reserve(mesh.triangles.size());
    for (const auto &mesh_tri: mesh.triangles) {
        FPolygon triangle(mesh_tri.nodes.size());
        for (int i = 0; i < mesh_tri.nodes.size(); ++i) {
            triangle[i] = mesh.point(mesh_tri.nodes[i]);
        }
        _triangles.push_back(std::move(triangle));
    }
    _mesh = std::move(mesh);
    _has_mesh = true;
    ComputeLimits(_triangles, _limits.x_min, _limits.x_max, _limits.y_min, _limits.y_max);
    if (map) {
        _map = std::move(*map);
        _has_map = true;
    } else {
        _has_map = false;
    }
    _has_pl = false;
}

void Trivis::FillPointLocationBuckets(
    std::optional<double> bucket_size,
    std::optional<double> max_avg_triangle_count_in_bucket
) {
    assert(_has_mesh);
    if (!max_avg_triangle_count_in_bucket) {
        _pl.Init(_limits, _triangles, bucket_size.value_or(ParamDefault::pl_bucket_size));
        _has_pl = true;
        return;
    }
    // Choose bucket size dynamically.
    auto avg_triangle_count = static_cast<double>(_triangles.size());
    double bucket_size_temp = std::max(_limits.x_max - _limits.x_min, _limits.y_max - _limits.y_min) / 2.0;
    while (avg_triangle_count > max_avg_triangle_count_in_bucket) {
        bucket_size_temp /= 2.0;
        _pl.Init(_limits, _triangles, bucket_size_temp);
        double triangles_per_bucket_total = 0.0;
        for (const auto &bucket: _pl.buckets().data()) {
            triangles_per_bucket_total += static_cast<double>(bucket.size());
        }
        avg_triangle_count = triangles_per_bucket_total / static_cast<double>(_pl.n_row() * _pl.n_col());
    }
    _has_pl = true;
}

void Trivis::OptimizePointLocationBucketTriangles() {
    assert(_has_pl);
    _pl.RemoveDuplicateTrianglesInBuckets();
    _pl.SortTrianglesInBucketsByLargestIntersectionArea(_triangles);
}

void Trivis::InitPointLocation(
    std::optional<double> bucket_size,
    std::optional<double> max_avg_triangle_count_in_bucket,
    std::optional<bool> optimize_bucket_triangles
) {
    assert(_has_mesh);
    FillPointLocationBuckets(bucket_size, max_avg_triangle_count_in_bucket);
    if (optimize_bucket_triangles.value_or(ParamDefault::pl_optimize_bucket_triangles)) {
        OptimizePointLocationBucketTriangles();
    }
}

void Trivis::SetPointLocationEpsilons(
    const std::optional<std::vector<double>> &eps1_seq,
    std::optional<double> eps2_squared
) {
    _pl_eps1_seq = eps1_seq.value_or(std::vector<double>{ParamDefault::pl_eps1.begin(), ParamDefault::pl_eps1.end()});
    _pl_eps2_squared = eps2_squared.value_or(ParamDefault::pl_eps2_squared);
}

void Trivis::Init(
    geom::PolyMap map,
    std::optional<double> pl_bucket_size,
    std::optional<double> pl_max_avg_triangle_count_in_bucket,
    std::optional<bool> pl_optimize_bucket_triangles,
    const std::optional<std::vector<double>> &pl_eps1_seq,
    std::optional<double> pl_eps2_squared
) {
    SetMap(std::move(map));
    ConstructMeshCDT();
    InitPointLocation(pl_bucket_size, pl_max_avg_triangle_count_in_bucket, pl_optimize_bucket_triangles.value_or(ParamDefault::pl_optimize_bucket_triangles));
    SetPointLocationEpsilons(pl_eps1_seq, pl_eps2_squared);
}

void Trivis::Init(
    mesh::TriMesh mesh, std::optional<geom::PolyMap> map,
    std::optional<double> pl_bucket_size,
    std::optional<double> pl_max_avg_triangle_count_in_bucket,
    std::optional<bool> pl_optimize_bucket_triangles,
    const std::optional<std::vector<double>> &pl_eps1_seq,
    std::optional<double> pl_eps2_squared
) {
    SetMesh(std::move(mesh), std::move(map));
    InitPointLocation(pl_bucket_size, pl_max_avg_triangle_count_in_bucket, pl_optimize_bucket_triangles.value_or(ParamDefault::pl_optimize_bucket_triangles));
    SetPointLocationEpsilons(pl_eps1_seq, pl_eps2_squared);
}

/// ##### LOCATE POINT ##### ///

std::optional<Trivis::PointLocationResult> Trivis::LocatePoint(
    const FPoint &q,
    const std::optional<std::vector<double>> &eps1_seq,
    const std::optional<double> &eps2_squared
) const {
    assert(_has_mesh);
    assert(_has_pl);
    const auto &eps1_seq_final = eps1_seq ? *eps1_seq : _pl_eps1_seq; // choose between given and saved options
    std::optional<int> tri_id_opt = _pl.FindTriangle(q, _triangles, eps1_seq_final);
    if (!tri_id_opt) {
        return std::nullopt;
    }
    double eps2_squared_final = eps2_squared.value_or(_pl_eps2_squared); // choose between given and saved options
    Trivis::PointLocationResult result;
    result.triangle_id = *tri_id_opt;
    const auto &q_triangle = _mesh.triangles[result.triangle_id];
    for (int node_id: q_triangle.nodes) {
        const auto &node_p = _mesh.point(node_id);
        if (q.SquaredDistanceTo(node_p) <= eps2_squared_final) {
            result.snap_to_nodes.push_back(node_id);
            if (q == node_p) {
                auto next_node = _mesh.nodes[node_id].next_weakly_intersect_node;
                if (next_node) {
                    while (next_node != node_id) {
                        result.snap_to_nodes.push_back(*next_node);
                        next_node = _mesh.nodes[*next_node].next_weakly_intersect_node;
                    }
                }
            }
            break;
        }
    }
    return result;
}

/// ##### VISIBILITY: INTERSECTION OF RAY AND OBSTACLE ##### ///

std::optional<Trivis::RayShootingResult> Trivis::ShootRay(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    const geom::FPoint &direction,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return ShootRay(q, q_location.triangle_id, direction, stats);
    } else {
        for (int node_id: q_location.snap_to_nodes) {
            auto ret = ShootRay(node_id, direction, stats);
            if (ret) {
                return ret;
            }
        }
        return std::nullopt;
    }
}

std::optional<Trivis::RayShootingResult> Trivis::ShootRay(
    const geom::FPoint &q,
    int q_triangle_id,
    const geom::FPoint &direction,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    const double w = _limits.x_max - _limits.x_min;
    const double h = _limits.y_max - _limits.y_min;
    const double multiplier = 2.0 * std::sqrt(w * w + h * h);
    const auto distant_t = q + (direction * multiplier);

    for (int edge_id: q_triangle.edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int node_l_id = edge.nodes[0];
        int node_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(node_r_id), _mesh.point(node_l_id))) {
            std::swap(node_l_id, node_r_id);
        }
        const auto &node_l_p = _mesh.point(node_l_id);
        const auto &node_r_p = _mesh.point(node_r_id);

        if (IsPointInCone(distant_t, node_l_p, q, node_r_p)) {
            if (edge.is_boundary()) {
                Trivis::RayShootingResult ret;
                ret.code = RaySegmentIntersection(q, distant_t, node_l_p, node_r_p, ret.p, ret.p2);
                if (ret.code == '0' || ret.code == 'c' || ret.code == 'e') {
                    return ret;
                }
                ret.edge_id = edge_id;
                if (ret.code == 'V' || ret.code == 's' || ret.code == 'v' || ret.code == 'i') {
                    if (ret.p == node_l_p) {
                        ret.node_id = node_l_id;
                        return ret;
                    }
                    if (ret.p == node_r_p) {
                        ret.node_id = node_r_id;
                        return ret;
                    }
                }
                return ret;
            }
            return ExpandEdgeObstacleIntersection(1, q, distant_t, node_l_id, node_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, stats);
        }
    }
    return std::nullopt;
}

std::optional<Trivis::RayShootingResult> Trivis::ShootRay(
    int node_id,
    const geom::FPoint &direction,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    if (node.edges.empty() || node.triangles.empty()) {
        return {};
    }

    const double w = _limits.x_max - _limits.x_min;
    const double h = _limits.y_max - _limits.y_min;
    const double multiplier = 2.0 * std::sqrt(w * w + h * h);
    const auto distant_t = node_p + (direction * multiplier);

    // for each triangle 'tri' adjacent to the node ...
    for (auto tri_id: node.triangles) {

        // for each edge of triangle 'tri' ...
        for (auto edge_id: _mesh.triangles[tri_id].edges) {

            const auto &edge = _mesh.edges[edge_id];
            int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

            int node_l_id = edge.nodes[0];
            int node_r_id = edge.nodes[1];

            // assume edge is opposite to node
            if (node_l_id != node_id && node_r_id != node_id) {

                // Establish the left and right restriction points.
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(node_r_id), _mesh.point(node_l_id))) {
                    std::swap(node_l_id, node_r_id);
                }
                const auto &node_l_p = _mesh.point(node_l_id);
                const auto &node_r_p = _mesh.point(node_r_id);

                if (IsPointInCone(distant_t, node_l_p, node_p, node_r_p)) {
                    if (edge.is_boundary()) {
                        Trivis::RayShootingResult ret;
                        ret.code = RaySegmentIntersection(node_p, distant_t, node_l_p, node_r_p, ret.p, ret.p2);
                        if (ret.code == '0' || ret.code == 'c' || ret.code == 'e') {
                            return ret;
                        }
                        ret.edge_id = edge_id;
                        if (ret.code == 'V' || ret.code == 's' || ret.code == 'v' || ret.code == 'i') {
                            if (ret.p == node_p) {
                                ret.node_id = node_id;
                                return ret;
                            }
                            if (ret.p == node_l_p) {
                                ret.node_id = node_l_id;
                                return ret;
                            }
                            if (ret.p == node_r_p) {
                                ret.node_id = node_r_id;
                                return ret;
                            }
                        }
                        return ret;
                    }
                    return ExpandEdgeObstacleIntersection(1, node_p, distant_t, node_l_id, node_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, stats);
                }
            }
        }
    }
    return std::nullopt;
}

/// ##### VISIBILITY: TWO-POINT QUERIES ##### ///

std::optional<bool> Trivis::IsVisible(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    const geom::FPoint &p,
    std::optional<double> radius,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return IsVisible(q, q_location.triangle_id, p, radius, stats);
    } else {
        bool at_least_one_value = false;
        for (int node_id: q_location.snap_to_nodes) {
            auto ret = IsVisible(node_id, p, radius, stats);
            if (ret.value_or(false)) {
                return true;
            }
            if (ret) {
                at_least_one_value = true;
            }
        }
        if (at_least_one_value) {
            return false;
        } else {
            return std::nullopt;
        }
    }
}

std::optional<bool> Trivis::IsVisible(
    const geom::FPoint &q,
    int q_triangle_id,
    const geom::FPoint &p,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    if (radius && q.SquaredDistanceTo(p) > (*radius) * (*radius)) {
        return false;
    }

    for (int edge_id: q_triangle.edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int node_l_id = edge.nodes[0];
        int node_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(node_r_id), _mesh.point(node_l_id))) {
            std::swap(node_l_id, node_r_id);
        }
        const auto &node_l_p = _mesh.point(node_l_id);
        const auto &node_r_p = _mesh.point(node_r_id);

        if (IsPointInCone(p, node_l_p, q, node_r_p)) {
            if (!TurnsLeft(node_l_p, node_r_p, p)) {
                // target lies in the same triangle as query
                return true;
            } else if (!edge.is_boundary()) {
                return ExpandEdgeVisibilityBetween(1, q, p, node_l_id, node_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, stats);
            }
        }
    }
    return std::nullopt;
}

std::optional<bool> Trivis::IsVisible(
    int node_id,
    const geom::FPoint &p,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    if (radius && node_p.SquaredDistanceTo(p) > (*radius) * (*radius)) {
        return false;
    }

    if (node.edges.empty() || node.triangles.empty()) {
        return false;
    }

    // for each triangle 'tri' adjacent to the node ...
    for (auto tri_id: node.triangles) {

        // for each edge of triangle 'tri' ...
        for (auto edge_id: _mesh.triangles[tri_id].edges) {

            const auto &edge = _mesh.edges[edge_id];
            int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

            int node_l_id = edge.nodes[0];
            int node_r_id = edge.nodes[1];

            // assume edge is opposite to node
            if (node_l_id != node_id && node_r_id != node_id) {

                // Establish the left and right restriction points.
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(node_r_id), _mesh.point(node_l_id))) {
                    std::swap(node_l_id, node_r_id);
                }
                const auto &node_l_p = _mesh.point(node_l_id);
                const auto &node_r_p = _mesh.point(node_r_id);

                if (IsPointInCone(p, node_l_p, node_p, node_r_p)) {
                    if (!TurnsLeft(node_l_p, node_r_p, p)) {
                        // target lies in the same triangle as query
                        return true;
                    } else if (!edge.is_boundary()) {
                        return ExpandEdgeVisibilityBetween(1, node_p, p, node_l_id, node_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, stats);
                    }
                }
            }
        }
    }
    return std::nullopt;
}

/// ##### VISIBILITY: MAP VERTICES ##### ///

std::optional<std::vector<int>> Trivis::VisibleVertices(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return VisibleVertices(q, q_location.triangle_id, nullptr, radius, stats);
    } else {
        std::optional<std::vector<int>> ret;
        for (int node_id: q_location.snap_to_nodes) {
            auto ret_opt = VisibleVertices(node_id, tabu_vertices, radius, stats);
            if (ret_opt) {
                if (!ret) {
                    ret = std::move(ret_opt);
                } else {
                    ret->insert(ret->end(), ret_opt->begin(), ret_opt->end());
                }
            }
        }
        return ret;
    }
}

std::optional<std::vector<int>> Trivis::VisibleVertices(
    const geom::FPoint &q,
    int q_triangle_id,
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    std::vector<int> ret;

    // Expand edges of the triangle.
    for (auto edge_id: _mesh.triangles[q_triangle_id].edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int rest_l_id = edge.nodes[0];
        int rest_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
            std::swap(rest_l_id, rest_r_id);
        }

        // Expand the edge.
        ExpandEdgeVisibleVertices(1, q, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, ret, tabu_vertices, stats);

    }

    return ret;
}

std::optional<std::vector<int>> Trivis::VisibleVertices(
    int node_id,
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    if (node.edges.empty() || node.triangles.empty()) {
        return {};
    }

    std::vector<int> ret;

    // Add the first edge.
    int e_id_1st = node.edges.front();
    const auto &e_1st = _mesh.edges[e_id_1st];
    int v_l_id_1st = e_1st.nodes[0];
    int v_r_id_1st = e_1st.nodes[1];
    if (_mesh.nodes[v_l_id_1st].edges.front() == e_id_1st) {
        std::swap(v_l_id_1st, v_r_id_1st);
    }
    if ((ret.empty() || ret.back() != v_r_id_1st) && (sq_radius < 0.0 || node_p.SquaredDistanceTo(_mesh.point(v_r_id_1st)) <= sq_radius)) {
        if (!tabu_vertices || !tabu_vertices->operator[](v_r_id_1st)) {
            ret.push_back(v_r_id_1st);
        }
    }
    if ((ret.empty() || ret.back() != v_l_id_1st) && (sq_radius < 0.0 || node_p.SquaredDistanceTo(_mesh.point(v_l_id_1st)) <= sq_radius)) {
        if (!tabu_vertices || !tabu_vertices->operator[](v_l_id_1st)) {
            ret.push_back(v_l_id_1st);
        }
    }

    for (int tri_id: node.triangles) {
        const auto &tri = _mesh.triangles[tri_id];
        for (int edge_id: tri.edges) {
            const auto &edge = _mesh.edges[edge_id];
            if (edge.nodes[0] != node_id && edge.nodes[1] != node_id) {
                int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

                // Establish the left and right restriction points.
                int rest_l_id = edge.nodes[0];
                int rest_r_id = edge.nodes[1];
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
                    std::swap(rest_l_id, rest_r_id);
                }

                // Expand the edge.
                ExpandEdgeVisibleVertices(1, node_p, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, ret, tabu_vertices, stats);

                break;
            }
        }

    }

    // Add the last edge.
    int e_id_lst = node.edges.back();
    const auto &e_lst = _mesh.edges[e_id_lst];
    int v_l_id_lst = e_lst.nodes[0];
    int v_r_id_lst = e_lst.nodes[1];
    if (_mesh.nodes[v_l_id_lst].edges.front() == e_id_lst) {
        std::swap(v_l_id_lst, v_r_id_lst);
    }
    if (v_r_id_lst != node_id && (ret.empty() || ret.back() != v_r_id_lst) && (sq_radius < 0.0 || node_p.SquaredDistanceTo(_mesh.point(v_r_id_lst)) <= sq_radius)) {
        if (!tabu_vertices || !tabu_vertices->operator[](v_r_id_lst)) {
            ret.push_back(v_r_id_lst);
        }
    }
    if (v_l_id_lst != node_id && (ret.empty() || ret.back() != v_l_id_lst) && (sq_radius < 0.0 || node_p.SquaredDistanceTo(_mesh.point(v_l_id_lst)) <= sq_radius)) {
        if (!tabu_vertices || !tabu_vertices->operator[](v_l_id_lst)) {
            ret.push_back(v_l_id_lst);
        }
    }

    return ret;
}

/// ##### VISIBILITY: INPUT POINTS ##### ///

std::optional<std::vector<int>> Trivis::VisiblePoints(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    std::optional<double> radius,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return VisiblePoints(q, q_location.triangle_id, points, points_locations, radius, stats);
    } else {
        std::optional<std::vector<int>> ret;
        for (int node_id: q_location.snap_to_nodes) {
            auto ret_opt = VisiblePoints(node_id, points, points_locations, radius, stats);
            if (ret_opt) {
                if (!ret) {
                    ret = std::move(ret_opt);
                } else {
                    ret->insert(ret->end(), ret_opt->begin(), ret_opt->end());
                }
            }
        }
        return ret;
    }
}

std::optional<std::vector<int>> Trivis::VisiblePoints(
    const geom::FPoint &q,
    int q_triangle_id,
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    std::vector<std::vector<int>> triangle_points(_mesh.triangles.size());
    for (int point_id = 0; point_id < points.size(); ++point_id) {
        const auto &point_location = points_locations[point_id];
        if (!point_location) {
            continue;
        }
        triangle_points[point_location->triangle_id].push_back(point_id);
    }

    return VisiblePoints(q, q_triangle_id, points, triangle_points, radius, stats);
}

std::optional<std::vector<int>> Trivis::VisiblePoints(
    int node_id,
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    std::vector<std::vector<int>> triangle_points(_mesh.triangles.size());
    for (int point_id = 0; point_id < points.size(); ++point_id) {
        const auto &point_location = points_locations[point_id];
        if (!point_location) {
            continue;
        }
        triangle_points[point_location->triangle_id].push_back(point_id);
    }

    return VisiblePoints(node_id, points, triangle_points, radius, stats);
}

std::optional<std::vector<int>> Trivis::VisiblePoints(
    const geom::FPoint &q,
    int q_triangle_id,
    const geom::FPoints &points,
    const std::vector<std::vector<int>> &triangle_points,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    std::vector<bool> point_visited(points.size(), false);

    std::vector<int> ret;

    // Add all points from the current q_triangle.
    for (int point_id: triangle_points[q_triangle_id]) {
        if (sq_radius < 0.0 || points[point_id].SquaredDistanceTo(q) <= sq_radius) {
            ret.push_back(point_id);
            point_visited[point_id] = true;
        }
    }

    // Expand edges of the q_triangle.
    for (auto edge_id: _mesh.triangles[q_triangle_id].edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int rest_l_id = edge.nodes[0];
        int rest_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
            std::swap(rest_l_id, rest_r_id);
        }

        // Expand the edge.
        ExpandEdgeVisiblePoints(points, triangle_points, 1, q, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, point_visited, ret, stats);

    }

    return ret;
}

std::optional<std::vector<int>> Trivis::VisiblePoints(
    int node_id,
    const geom::FPoints &points,
    const std::vector<std::vector<int>> &triangle_points,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    if (node.edges.empty() || node.triangles.empty()) {
        return {};
    }

    std::vector<bool> point_visited(points.size(), false);

    std::vector<int> ret;

    for (int tri_id: node.triangles) {
        const auto &tri = _mesh.triangles[tri_id];

        // Add all points from the current triangle.
        for (int point_id: triangle_points[tri_id]) {
            if (sq_radius < 0.0 || points[point_id].SquaredDistanceTo(node_p) <= sq_radius) {
                ret.push_back(point_id);
                point_visited[point_id] = true;
            }
        }

        for (int edge_id: tri.edges) {
            const auto &edge = _mesh.edges[edge_id];
            if (edge.nodes[0] != node_id && edge.nodes[1] != node_id) {
                int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

                // Establish the left and right restriction points.
                int rest_l_id = edge.nodes[0];
                int rest_r_id = edge.nodes[1];
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
                    std::swap(rest_l_id, rest_r_id);
                }

                // Expand the edge.
                ExpandEdgeVisiblePoints(points, triangle_points, 1, node_p, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, point_visited, ret, stats);

                break;
            }
        }

    }

    return ret;
}

/// ##### VISIBILITY: VISIBILITY GRAPHS ##### ///

std::optional<std::vector<std::vector<int>>> Trivis::VertexVertexVisibilityGraph(
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }
    ExpansionStats stats_temp;
    ExpansionStats *stats_temp_ptr = stats ? &stats_temp : nullptr;

    std::vector<std::vector<int>> ret(_mesh.nodes.size());
    for (int node_id = 0; node_id < _mesh.nodes.size(); ++node_id) {
        if (tabu_vertices && tabu_vertices->operator[](node_id)) {
            continue;
        }
        auto visible_vertices_opt = VisibleVertices(node_id, tabu_vertices, radius, stats_temp_ptr);
        if (!visible_vertices_opt) {
            continue;
        }
        ret[node_id] = std::move(*visible_vertices_opt);
        if (stats) {
            stats->num_expansions += stats_temp_ptr->num_expansions;
            stats->max_recursion_depth = std::max(stats->max_recursion_depth, stats_temp_ptr->max_recursion_depth);
        }
    }
    return ret;
}

std::optional<std::vector<std::vector<bool>>> Trivis::VertexVertexVisibilityGraphBool(
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }
    ExpansionStats stats_temp;
    ExpansionStats *stats_temp_ptr = stats ? &stats_temp : nullptr;

    std::vector<std::vector<bool>> ret(_mesh.nodes.size(), std::vector<bool>(_mesh.nodes.size(), false));
    for (int node_id = 0; node_id < _mesh.nodes.size(); ++node_id) {
        ret[node_id][node_id] = true;
        if (tabu_vertices && tabu_vertices->operator[](node_id)) {
            continue;
        }
        auto visible_vertices_opt = VisibleVertices(node_id, tabu_vertices, radius, stats_temp_ptr);
        if (!visible_vertices_opt) {
            continue;
        }
        for (int visible_node_id: *visible_vertices_opt) {
            ret[node_id][visible_node_id] = ret[visible_node_id][node_id] = true;
        }
        if (stats) {
            stats->num_expansions += stats_temp_ptr->num_expansions;
            stats->max_recursion_depth = std::max(stats->max_recursion_depth, stats_temp_ptr->max_recursion_depth);
        }
    }
    return ret;
}

std::optional<std::vector<std::vector<int>>> Trivis::VertexPointVisibilityGraph(
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }
    ExpansionStats stats_temp;
    ExpansionStats *stats_temp_ptr = stats ? &stats_temp : nullptr;

    std::vector<std::vector<int>> ret(points.size());
    for (int point_id = 0; point_id < points.size(); ++point_id) {
        const auto &p = points[point_id];
        const auto &p_location = points_locations[point_id];
        if (!p_location) {
            continue;
        }
        std::vector<int> visible_vertices;
        if (!p_location->snap_to_nodes.empty()) {
            for (int node_id: p_location->snap_to_nodes) {
                auto visible_vertices_opt = VisibleVertices(node_id, tabu_vertices, radius, stats);
                if (!visible_vertices_opt) {
                    continue;
                }
                visible_vertices.insert(visible_vertices.begin(), visible_vertices_opt->begin(), visible_vertices_opt->end());
            }
        } else {
            auto visible_vertices_opt = VisibleVertices(p, p_location->triangle_id, tabu_vertices, radius, stats);
            if (!visible_vertices_opt) {
                continue;
            }
            visible_vertices = std::move(*visible_vertices_opt);
        }
        ret[point_id] = std::move(visible_vertices);
        if (stats) {
            stats->num_expansions += stats_temp_ptr->num_expansions;
            stats->max_recursion_depth = std::max(stats->max_recursion_depth, stats_temp_ptr->max_recursion_depth);
        }
    }

    return ret;
}

std::optional<std::vector<std::vector<bool>>> Trivis::VertexPointVisibilityGraphBool(
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    const std::vector<bool> *tabu_vertices,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }
    ExpansionStats stats_temp;
    ExpansionStats *stats_temp_ptr = stats ? &stats_temp : nullptr;

    std::vector<std::vector<bool>> ret(points.size(), std::vector<bool>(_mesh.nodes.size(), false));
    for (int point_id = 0; point_id < points.size(); ++point_id) {
        const auto &p = points[point_id];
        const auto &p_location = points_locations[point_id];
        if (!p_location) {
            continue;
        }
        std::vector<int> visible_vertices;
        if (!p_location->snap_to_nodes.empty()) {
            for (int node_id: p_location->snap_to_nodes) {
                auto visible_vertices_opt = VisibleVertices(node_id, tabu_vertices, radius, stats);
                if (!visible_vertices_opt) {
                    continue;
                }
                visible_vertices.insert(visible_vertices.begin(), visible_vertices_opt->begin(), visible_vertices_opt->end());
            }
        } else {
            auto visible_vertices_opt = VisibleVertices(p, p_location->triangle_id, tabu_vertices, radius, stats);
            if (!visible_vertices_opt) {
                continue;
            }
            visible_vertices = std::move(*visible_vertices_opt);
        }
        for (int node_id: visible_vertices) {
            ret[point_id][node_id] = true;
        }
        if (stats) {
            stats->num_expansions += stats_temp_ptr->num_expansions;
            stats->max_recursion_depth = std::max(stats->max_recursion_depth, stats_temp_ptr->max_recursion_depth);
        }
    }

    return ret;
}

std::optional<std::vector<std::vector<int>>> Trivis::PointPointVisibilityGraph(
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }
    ExpansionStats stats_temp;
    ExpansionStats *stats_temp_ptr = stats ? &stats_temp : nullptr;

    std::vector<std::vector<int>> ret(points.size());
    for (int point_id = 0; point_id < points.size(); ++point_id) {
        const auto &p = points[point_id];
        const auto &p_location = points_locations[point_id];
        if (!p_location) {
            continue;
        }
        std::vector<int> visible_points;
        if (!p_location->snap_to_nodes.empty()) {
            for (int node_id: p_location->snap_to_nodes) {
                auto visible_points_opt = VisiblePoints(node_id, points, points_locations, radius, stats);
                if (!visible_points_opt) {
                    continue;
                }
                visible_points.insert(visible_points.begin(), visible_points_opt->begin(), visible_points_opt->end());
            }
        } else {
            auto visible_points_opt = VisiblePoints(p, p_location->triangle_id, points, points_locations, radius, stats);
            if (!visible_points_opt) {
                continue;
            }
            visible_points = std::move(*visible_points_opt);
        }
        ret[point_id] = std::move(visible_points);
        if (stats) {
            stats->num_expansions += stats_temp_ptr->num_expansions;
            stats->max_recursion_depth = std::max(stats->max_recursion_depth, stats_temp_ptr->max_recursion_depth);
        }
    }

    return ret;
}

std::optional<std::vector<std::vector<bool>>> Trivis::PointPointVisibilityGraphBool(
    const geom::FPoints &points,
    const std::vector<std::optional<PointLocationResult>> &points_locations,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }
    ExpansionStats stats_temp;
    ExpansionStats *stats_temp_ptr = stats ? &stats_temp : nullptr;

    std::vector<std::vector<bool>> ret(points.size(), std::vector<bool>(points.size(), false));
    for (int point_id = 0; point_id < points.size(); ++point_id) {
        ret[point_id][point_id] = true;
        const auto &p = points[point_id];
        const auto &p_location = points_locations[point_id];
        if (!p_location) {
            continue;
        }
        std::vector<int> visible_points;
        if (!p_location->snap_to_nodes.empty()) {
            for (int node_id: p_location->snap_to_nodes) {
                auto visible_points_opt = VisiblePoints(node_id, points, points_locations, radius, stats);
                if (!visible_points_opt) {
                    continue;
                }
                visible_points.insert(visible_points.begin(), visible_points_opt->begin(), visible_points_opt->end());
            }
        } else {
            auto visible_points_opt = VisiblePoints(p, p_location->triangle_id, points, points_locations, radius, stats);
            if (!visible_points_opt) {
                continue;
            }
            visible_points = std::move(*visible_points_opt);
        }
        for (int visible_point_id: visible_points) {
            ret[point_id][visible_point_id] = ret[visible_point_id][point_id] = true;
        }
        if (stats) {
            stats->num_expansions += stats_temp_ptr->num_expansions;
            stats->max_recursion_depth = std::max(stats->max_recursion_depth, stats_temp_ptr->max_recursion_depth);
        }
    }

    return ret;
}


/// ##### VISIBILITY: VISIBILITY REGIONS (ABSTRACT REPRESENTATION) ##### ///

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegion(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    std::optional<double> radius,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return VisibilityRegion(q, q_location.triangle_id, radius, stats);
    } else {
        std::optional<AbstractVisibilityRegion> ret;
        for (int node_id: q_location.snap_to_nodes) {
            auto ret_opt = VisibilityRegion(node_id, radius, stats);
            if (ret_opt) {
                if (!ret) {
                    ret = std::move(ret_opt);
                } else {
                    ret->segments.insert(ret->segments.end(), ret_opt->segments.begin(), ret_opt->segments.end());
                }
            }
        }
        return ret;
    }
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegion(
    const geom::FPoint &q,
    int q_triangle_id,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Use the node version in case q is a node of the triangle.
    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    // Save the seed.
    AbstractVisibilityRegion ret;
    ret.seed_id = -1; // seed is NOT one of mesh's nodes
    ret.seed = q;
    ret.segments.reserve(mesh().nodes.size());

    // Expand edges of the triangle.
    for (auto edge_id: _mesh.triangles[q_triangle_id].edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int rest_l_id = edge.nodes[0];
        int rest_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
            std::swap(rest_l_id, rest_r_id);
        }

        // Expand the edge.
        ExpandEdgeVisibilityRegion(1, q, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, ret, stats);

    }

    return ret;
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegion(
    int node_id,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    // Save the seed.
    AbstractVisibilityRegion ret;
    ret.seed_id = node_id; // seed is one of mesh's nodes
    ret.seed = node_p;

    if (node.edges.empty() || node.triangles.empty()) {
        return ret;
    }

    ret.segments.reserve(mesh().nodes.size());

    // Add the first edge.
    int e_id_1st = node.edges.front();
    const auto &e_1st = _mesh.edges[e_id_1st];
    AbstractVisibilityRegionSegment seg_1st;
    int v_l_id_1st = e_1st.nodes[0];
    int v_r_id_1st = e_1st.nodes[1];
    if (_mesh.nodes[v_l_id_1st].edges.front() == e_id_1st) {
        std::swap(v_l_id_1st, v_r_id_1st);
    }
    seg_1st.id = e_id_1st;
    seg_1st.v1.id = v_r_id_1st;
    seg_1st.v2.id = v_l_id_1st;

    ret.segments.push_back(seg_1st);

    for (int tri_id: node.triangles) {
        const auto &tri = _mesh.triangles[tri_id];
        for (int edge_id: tri.edges) {

            const auto &edge = _mesh.edges[edge_id];
            if (edge.nodes[0] != node_id && edge.nodes[1] != node_id) {
                int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

                // Establish the left and right restriction points.
                int rest_l_id = edge.nodes[0];
                int rest_r_id = edge.nodes[1];
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
                    std::swap(rest_l_id, rest_r_id);
                }

                // Expand the edge.
                ExpandEdgeVisibilityRegion(1, node_p, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, ret, stats);

                break;
            }
        }

    }

    // Add the last edge.
    int e_id_lst = node.edges.back();
    const auto &e_lst = _mesh.edges[e_id_lst];
    AbstractVisibilityRegionSegment seg_lst;
    int v_l_id_lst = e_lst.nodes[0];
    int v_r_id_lst = e_lst.nodes[1];
    if (_mesh.nodes[v_l_id_lst].edges.front() == e_id_lst) {
        std::swap(v_l_id_lst, v_r_id_lst);
    }
    seg_lst.id = e_id_lst;
    seg_lst.v1.id = v_r_id_lst;
    seg_lst.v2.id = v_l_id_lst;

    ret.segments.push_back(seg_lst);

    return ret;
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegionIterative(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    std::optional<double> radius,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return VisibilityRegionIterative(q, q_location.triangle_id, radius, stats);
    } else {
        std::optional<AbstractVisibilityRegion> ret;
        for (int node_id: q_location.snap_to_nodes) {
            auto ret_opt = VisibilityRegionIterative(node_id, radius, stats);
            if (ret_opt) {
                if (!ret) {
                    ret = std::move(ret_opt);
                } else {
                    ret->segments.insert(ret->segments.end(), ret_opt->segments.begin(), ret_opt->segments.end());
                }
            }
        }
        return ret;
    }
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegionIterative(
    const geom::FPoint &q,
    int q_triangle_id,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Use the node version in case q is a node of the triangle.
    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    // Save the seed.
    AbstractVisibilityRegion ret;
    ret.seed_id = -1; // seed is NOT one of mesh's nodes
    ret.seed = q;
    ret.segments.reserve(mesh().nodes.size());

    bool has_expansion1 = false;
    Trivis::EdgeExpansionInfo expansion1;
    bool has_expansion2 = false;
    Trivis::EdgeExpansionInfo expansion2;
    std::vector<Trivis::EdgeExpansionInfo> stack;
    bool next_expansion_valid = false;
    Trivis::EdgeExpansionInfo next_expansion;


    // Expand edges of the triangle.
    for (auto edge_id: q_triangle.edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int rest_l_id = edge.nodes[0];
        int rest_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
            std::swap(rest_l_id, rest_r_id);
        }

        next_expansion.level = 1;
        next_expansion.rest_l_id = rest_l_id;
        next_expansion.rest_r_id = rest_r_id;
        next_expansion.curr_edge_id = edge_id;
        next_expansion.curr_edge_tri_id = edge_tri_id == 0 ? 1 : 0;
        next_expansion_valid = true;

        // Iterative depth-first search instead of recursion.
        while (true) {
            if (next_expansion_valid) {
                ExpandEdgeVisibilityRegionIterative(q, next_expansion, sq_radius, ret, has_expansion1, expansion1, has_expansion2, expansion2, stats);
            } else {
                ExpandEdgeVisibilityRegionIterative(q, stack.back(), sq_radius, ret, has_expansion1, expansion1, has_expansion2, expansion2, stats);
                stack.pop_back();
            }
            if (has_expansion1) {
                if (has_expansion2) {
                    stack.push_back(expansion2);
                }
                next_expansion = expansion1;
                next_expansion_valid = true;
                continue;
            }
            if (has_expansion2) {
                next_expansion = expansion2;
                next_expansion_valid = true;
                continue;
            }
            if (stack.empty()) {
                break;
            }
            next_expansion_valid = false;
        }
    }

    return ret;
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegionIterative(
    int node_id,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    // Save the seed.
    AbstractVisibilityRegion ret;
    ret.seed_id = -1; // seed is NOT one of mesh's nodes
    ret.seed = node_p;

    if (node.edges.empty() || node.triangles.empty()) {
        return ret;
    }

    ret.segments.reserve(mesh().nodes.size());

    // Add the first edge.
    int e_id_1st = node.edges.front();
    const auto &e_1st = _mesh.edges[e_id_1st];
    AbstractVisibilityRegionSegment seg_1st;
    int v_l_id_1st = e_1st.nodes[0];
    int v_r_id_1st = e_1st.nodes[1];
    if (_mesh.nodes[v_l_id_1st].edges.front() == e_id_1st) {
        std::swap(v_l_id_1st, v_r_id_1st);
    }
    seg_1st.id = e_id_1st;
    seg_1st.v1.id = v_r_id_1st;
    seg_1st.v2.id = v_l_id_1st;

    ret.segments.push_back(seg_1st);

    bool has_expansion1 = false;
    Trivis::EdgeExpansionInfo expansion1;
    bool has_expansion2 = false;
    Trivis::EdgeExpansionInfo expansion2;
    std::vector<Trivis::EdgeExpansionInfo> stack;

    for (int tri_id: node.triangles) {
        const auto &tri = _mesh.triangles[tri_id];
        for (int edge_id: tri.edges) {

            const auto &edge = _mesh.edges[edge_id];
            if (edge.nodes[0] != node_id && edge.nodes[1] != node_id) {
                int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

                // Establish the left and right restriction points.
                int rest_l_id = edge.nodes[0];
                int rest_r_id = edge.nodes[1];
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
                    std::swap(rest_l_id, rest_r_id);
                }

                // Iterative depth-first search instead of recursion.
                stack.push_back({1, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0});
                while (!stack.empty()) {
                    ExpandEdgeVisibilityRegionIterative(node_p, stack.back(), sq_radius, ret, has_expansion1, expansion1, has_expansion2, expansion2, stats);
                    stack.pop_back();
                    if (has_expansion2) {
                        stack.push_back(expansion2);
                    }
                    if (has_expansion1) {
                        stack.push_back(expansion1);
                    }
                }

                break;
            }
        }

    }

    // Add the last edge.
    int e_id_lst = node.edges.back();
    const auto &e_lst = _mesh.edges[e_id_lst];
    AbstractVisibilityRegionSegment seg_lst;
    int v_l_id_lst = e_lst.nodes[0];
    int v_r_id_lst = e_lst.nodes[1];
    if (_mesh.nodes[v_l_id_lst].edges.front() == e_id_lst) {
        std::swap(v_l_id_lst, v_r_id_lst);
    }
    seg_lst.id = e_id_lst;
    seg_lst.v1.id = v_r_id_lst;
    seg_lst.v2.id = v_l_id_lst;

    ret.segments.push_back(seg_lst);

    return ret;
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegionWithHistory(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    std::vector<ExpansionHistoryStep> &history,
    std::optional<double> radius,
    ExpansionStats *stats
) const {
    if (q_location.snap_to_nodes.empty()) {
        return VisibilityRegionWithHistory(q, q_location.triangle_id, history, radius, stats);
    } else {
        std::optional<AbstractVisibilityRegion> ret;
        for (int node_id: q_location.snap_to_nodes) {
            auto ret_opt = VisibilityRegionWithHistory(node_id, history, radius, stats);
            if (ret_opt) {
                if (!ret) {
                    ret = std::move(ret_opt);
                } else {
                    ret->segments.insert(ret->segments.end(), ret_opt->segments.begin(), ret_opt->segments.end());
                }
            }
        }
        return ret;
    }
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegionWithHistory(
    const geom::FPoint &q,
    int q_triangle_id,
    std::vector<ExpansionHistoryStep> &history,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(0 <= q_triangle_id && q_triangle_id < _mesh.triangles.size());
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    const auto &q_triangle = _mesh.triangles[q_triangle_id];

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    // Save the seed.
    AbstractVisibilityRegion ret;
    ret.seed_id = -1; // seed is NOT one of mesh's nodes
    ret.seed = q;
    ret.segments.reserve(mesh().nodes.size());

    // Expand edges of the triangle.
    for (auto edge_id: _mesh.triangles[q_triangle_id].edges) {
        const auto &edge = _mesh.edges[edge_id];
        int edge_tri_id = edge.triangles[0] == q_triangle_id ? 0 : 1;

        // Establish the left and right restriction points.
        int rest_l_id = edge.nodes[0];
        int rest_r_id = edge.nodes[1];
        if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
            std::swap(rest_l_id, rest_r_id);
        }

        // Expand the edge.
        ExpandEdgeVisibilityRegionWithHistory(1, q, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, ret, history, stats);

    }

    return ret;
}

std::optional<AbstractVisibilityRegion> Trivis::VisibilityRegionWithHistory(
    int node_id,
    std::vector<ExpansionHistoryStep> &history,
    std::optional<double> radius,
    ExpansionStats *stats
) const {

    assert(_has_mesh);
    assert(!radius || radius > 0.0);

    if (stats) {
        stats->num_expansions = 0;
        stats->max_recursion_depth = 0;
    }

    // Compute radius squared (or set to -1 in case it is not given).
    double sq_radius = (radius && radius > 0.0) ? *radius * *radius : -1.0;

    const auto &node = _mesh.nodes[node_id];
    const auto &node_p = node.point;

    // Save the seed.
    AbstractVisibilityRegion ret;
    ret.seed_id = node_id; // seed is one of mesh's nodes
    ret.seed = node_p;

    if (node.edges.empty() || node.triangles.empty()) {
        return ret;
    }

    ret.segments.reserve(mesh().nodes.size());

    // Add the first edge.
    int e_id_1st = node.edges.front();
    const auto &e_1st = _mesh.edges[e_id_1st];
    AbstractVisibilityRegionSegment seg_1st;
    int v_l_id_1st = e_1st.nodes[0];
    int v_r_id_1st = e_1st.nodes[1];
    if (_mesh.nodes[v_l_id_1st].edges.front() == e_id_1st) {
        std::swap(v_l_id_1st, v_r_id_1st);
    }
    seg_1st.id = e_id_1st;
    seg_1st.v1.id = v_r_id_1st;
    seg_1st.v2.id = v_l_id_1st;

    ret.segments.push_back(seg_1st);

    for (int tri_id: node.triangles) {
        const auto &tri = _mesh.triangles[tri_id];
        for (int edge_id: tri.edges) {

            const auto &edge = _mesh.edges[edge_id];
            if (edge.nodes[0] != node_id && edge.nodes[1] != node_id) {
                int edge_tri_id = edge.triangles[0] == tri_id ? 0 : 1;

                // Establish the left and right restriction points.
                int rest_l_id = edge.nodes[0];
                int rest_r_id = edge.nodes[1];
                if (!TurnsLeft(_mesh.point(edge.opposites[edge_tri_id]), _mesh.point(rest_r_id), _mesh.point(rest_l_id))) {
                    std::swap(rest_l_id, rest_r_id);
                }

                // Expand the edge.
                ExpandEdgeVisibilityRegionWithHistory(1, node_p, rest_l_id, rest_r_id, edge_id, edge_tri_id == 0 ? 1 : 0, sq_radius, ret, history, stats);

                break;
            }
        }

    }

    // Add the last edge.
    int e_id_lst = node.edges.back();
    const auto &e_lst = _mesh.edges[e_id_lst];
    AbstractVisibilityRegionSegment seg_lst;
    int v_l_id_lst = e_lst.nodes[0];
    int v_r_id_lst = e_lst.nodes[1];
    if (_mesh.nodes[v_l_id_lst].edges.front() == e_id_lst) {
        std::swap(v_l_id_lst, v_r_id_lst);
    }
    seg_lst.id = e_id_lst;
    seg_lst.v1.id = v_r_id_lst;
    seg_lst.v2.id = v_l_id_lst;

    ret.segments.push_back(seg_lst);

    return ret;
}

/// ##### VISIBILITY REGIONS POSTPROCESSING AND UTILITIES ##### ///

RadialVisibilityRegion Trivis::ComputeIntersections(
    const AbstractVisibilityRegion &abstract,
    bool line_line_mode
) const {
    assert(_has_mesh);
    RadialVisibilityRegion ret;
    ret.seed_id = abstract.seed_id;
    ret.seed = abstract.seed;
    int n_minus_1 = static_cast<int>(abstract.segments.size()) - 1;
    for (int i_prev = n_minus_1, i = 0; i < abstract.segments.size(); i_prev = i++) {
        AbstractVisibilityRegionSegment seg_prev = abstract.segments[i_prev];
        AbstractVisibilityRegionSegment seg = abstract.segments[i];
        if (seg_prev.v2.is_intersection || seg_prev.v2.id != seg.v1.id) {
            // Append v1.
            AppendNotCollinearWithIntersection(abstract.seed, seg.v1, -1, false, line_line_mode, ret);
        }
        // Append v2.
        AppendNotCollinearWithIntersection(abstract.seed, seg.v2, seg.id, i == n_minus_1, line_line_mode, ret);
    }
    return ret;
}

bool Trivis::IsValid(const RadialVisibilityRegion &visibility_region) {
    return IsValid(visibility_region);
}

void Trivis::RemoveAntennas(
    RadialVisibilityRegion &visibility_region
) {
    assert(IsValid(visibility_region));
    RemoveAntennas(visibility_region);
}

RadialVisibilityRegion Trivis::IntersectWithCircle(
    double radius,
    const RadialVisibilityRegion &visibility_region
) {
    assert(IsValid(visibility_region));
    RadialVisibilityRegion ret;
    ret.radius = radius;
    ret.seed_id = visibility_region.seed_id;
    ret.seed = visibility_region.seed;
    double sq_radius = radius * radius;
    int n_minus_1 = static_cast<int>(visibility_region.vertices.size()) - 1;
    for (int i_prev = n_minus_1, i = 0; i < visibility_region.vertices.size(); i_prev = i++) {
        const auto &vi_prev = visibility_region.vertices[i_prev];
        const auto &vi = visibility_region.vertices[i];
        if (vi.edge_flag < -1) {
            // The whole edge is outside the radius: IGNORE IT.
            continue;
        }
        bool is_inside_pi_prev = vi_prev.point.SquaredDistanceTo(visibility_region.seed) <= sq_radius;
        bool is_inside_pi = vi.point.SquaredDistanceTo(visibility_region.seed) <= sq_radius;
        if (is_inside_pi_prev) {
            if (is_inside_pi) {
                // The whole edge is inside the radius (no intersection).
                if (vi_prev.edge_flag < -1) {
                    ret.vertices.push_back({vi_prev.vertex_flag, -2, vi_prev.point});
                }
                ret.vertices.push_back({vi.vertex_flag, vi.edge_flag, vi.point});
            } else {
                // The first endpoint is inside, the second is outside (1 intersection).
                FPoint intersection;
                if (vi.edge_flag < 0) {
                    intersection = vi_prev.point - visibility_region.seed;
                    intersection = intersection / intersection.Norm();
                    intersection = visibility_region.seed + intersection * radius;
                } else {
                    auto intersections = LineCircleIntersections(vi_prev.point, vi.point, visibility_region.seed, radius, true);
                    if (intersections.empty()) {
                        intersections = LineCircleIntersections(vi_prev.point.CopySwappedXY(), vi.point.CopySwappedXY(), visibility_region.seed.CopySwappedXY(), radius, true);
                        if (intersections.empty()) {
                            std::cerr << std::fixed << "[TriVis][IntersectWithCircle] LineCircle intersections should not be empty (1)! ";
                            std::cerr << "Line: " << vi_prev.point << ", " << vi.point << ", Circle:" << visibility_region.seed << ", " << radius << ".\n";
                            continue;
                        }
                        intersections[0].SwapXY();
                    }
                    intersection = intersections[0];
                }
                ret.vertices.push_back({-1, vi.edge_flag, intersection});
            }
        } else if (is_inside_pi) {
            // The first endpoint is outside, the second is inside (1 intersection).
            FPoint intersection;
            if (vi.edge_flag < 0) {
                intersection = vi_prev.point - visibility_region.seed;
                intersection = intersection / intersection.Norm();
                intersection = visibility_region.seed + intersection * radius;
            } else {
                auto intersections = LineCircleIntersections(vi_prev.point, vi.point, visibility_region.seed, radius, true);
                if (intersections.empty()) {
                    intersections = LineCircleIntersections(vi_prev.point.CopySwappedXY(), vi.point.CopySwappedXY(), visibility_region.seed.CopySwappedXY(), radius, true);
                    if (intersections.empty()) {
                        std::cerr << std::fixed << "[TriVis][IntersectWithCircle] LineCircle intersections should not be empty (2)! ";
                        std::cerr << "Line: " << vi_prev.point << ", " << vi.point << ", Circle:" << visibility_region.seed << ", " << radius << ".\n";
                        continue;
                    }
                    intersections[0].SwapXY();
                }
                intersection = intersections[0];
            }
            ret.vertices.push_back({-1, -2, intersection});
            ret.vertices.push_back(vi);
        } else {
            // Both endpoints are outside (2 intersections).
            auto intersections = LineCircleIntersections(vi_prev.point, vi.point, visibility_region.seed, radius, true);
            if (!intersections.empty()) {
                if (intersections.size() > 1 && vi_prev.point.SquaredDistanceTo(intersections[0]) > vi_prev.point.SquaredDistanceTo(intersections[1])) {
                    std::swap(intersections[0], intersections[1]);
                }
                ret.vertices.push_back({-1, -2, intersections[0]});
                if (intersections.size() > 1) {
                    ret.vertices.push_back({-1, vi.edge_flag, intersections[1]});
                }
            }
        }
    }
    return ret;
}

void Trivis::RemoveShortEdges(
    double min_edge_length,
    RadialVisibilityRegion &visibility_region
) {
    assert(IsValid(visibility_region));
    RemoveShortEdges(min_edge_length, visibility_region);
}

RadialVisibilityRegion Trivis::SampleArcEdges(
    const RadialVisibilityRegion &visibility_region,
    double max_angle
) {
    assert(IsValid(visibility_region));
    return SampleArcEdges(visibility_region, max_angle);
}

RadialVisibilityRegion Trivis::Postprocess(
    const AbstractVisibilityRegion &abstract,
    bool line_line_mode_intersections,
    bool remove_antennas,
    std::optional<double> radius,
    std::optional<double> min_edge_length,
    std::optional<double> sampling_max_angle
) const {
    auto ret = ComputeIntersections(abstract, line_line_mode_intersections);
    if (remove_antennas) {
        RemoveAntennas(ret);
    }
    if (radius) {
        ret = IntersectWithCircle(*radius, ret);
    }
    if (min_edge_length) {
        RemoveShortEdges(*min_edge_length, ret);
    }
    if (sampling_max_angle) {
        ret = SampleArcEdges(ret, *sampling_max_angle);
    }
    return ret;
}

geom::FPolygon Trivis::ConvertToPolygon(
    const RadialVisibilityRegion &visibility_region
) {
    assert(IsValid(visibility_region));
    geom::FPolygon ret;
    ret.reserve(visibility_region.vertices.size());
    for (const auto &v: visibility_region.vertices) {
        ret.push_back(v.point);
    }
    return ret;
}

/// ##### SHORTCUTS: VISIBILITY REGIONS + POSTPROCESSING ##### ///

std::optional<RadialVisibilityRegion> Trivis::VisibilityRegionWithPostprocessing(
    const geom::FPoint &q,
    const PointLocationResult &q_location,
    std::optional<double> radius,
    bool line_line_mode_intersections,
    bool remove_antennas,
    std::optional<double> min_edge_length,
    std::optional<double> sampling_max_angle,
    ExpansionStats *stats
) const {
    auto abstract_opt = VisibilityRegion(q, q_location, radius, stats);
    if (!abstract_opt) {
        return std::nullopt;
    }
    return Postprocess(*abstract_opt, line_line_mode_intersections, remove_antennas, radius, min_edge_length, sampling_max_angle);
}

std::optional<RadialVisibilityRegion> Trivis::VisibilityRegionWithPostprocessing(
    const geom::FPoint &q,
    int q_triangle_id,
    std::optional<double> radius,
    bool line_line_mode_intersections,
    bool remove_antennas,
    std::optional<double> min_edge_length,
    std::optional<double> sampling_max_angle,
    ExpansionStats *stats
) const {
    auto abstract_opt = VisibilityRegion(q, q_triangle_id, radius, stats);
    if (!abstract_opt) {
        return std::nullopt;
    }
    return Postprocess(*abstract_opt, line_line_mode_intersections, remove_antennas, radius, min_edge_length, sampling_max_angle);
}

std::optional<RadialVisibilityRegion> Trivis::VisibilityRegionWithPostprocessing(
    int node_id,
    std::optional<double> radius,
    bool line_line_mode_intersections,
    bool remove_antennas,
    std::optional<double> min_edge_length,
    std::optional<double> sampling_max_angle,
    ExpansionStats *stats
) const {
    auto abstract_opt = VisibilityRegion(node_id, radius, stats);
    if (!abstract_opt) {
        return std::nullopt;
    }
    return Postprocess(*abstract_opt, line_line_mode_intersections, remove_antennas, radius, min_edge_length, sampling_max_angle);
}

/// ##### GENERAL UTILITIES ##### ///

bool Trivis::SaveMesh(
    const mesh::TriMesh &mesh,
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
    fs << "\n[NODES]\n";
    fs << mesh.nodes.size() << "\n";
    for (int node_id = 0; node_id < mesh.nodes.size(); ++node_id) {
        const auto &node = mesh.nodes[node_id];
        fs << node_id;
        // point
        fs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << node.point.x;
        fs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << node.point.y;
        // edges
        fs << " " << node.edges.size();
        for (int edge_id: node.edges) {
            fs << " " << edge_id;
        }
        // triangles
        fs << " " << node.triangles.size();
        for (int triangle_id: node.triangles) {
            fs << " " << triangle_id;
        }
        // split partner
        fs << " " << node.next_weakly_intersect_node.value_or(-1);
        fs << "\n";
    }
    fs << "\n[EDGES]\n";
    fs << mesh.edges.size() << "\n";
    for (int edge_id = 0; edge_id < mesh.edges.size(); ++edge_id) {
        const auto &edge = mesh.edges[edge_id];
        fs << edge_id;
        // nodes
        fs << " " << edge.nodes[0] << " " << edge.nodes[1];
        // triangles
        fs << " " << edge.triangles.size();
        for (int triangle_id: edge.triangles) {
            fs << " " << triangle_id;
        }
        // opposites
        fs << " " << edge.opposites.size();
        for (int node_id: edge.opposites) {
            fs << " " << node_id;
        }
        fs << "\n";
    }

    fs << "\n[TRIANGLES]\n";
    fs << mesh.triangles.size() << "\n";
    for (int triangle_id = 0; triangle_id < mesh.triangles.size(); ++triangle_id) {
        const auto &triangle = mesh.triangles[triangle_id];
        fs << triangle_id;
        // nodes
        fs << " " << triangle.nodes[0] << " " << triangle.nodes[1] << " " << triangle.nodes[2];
        // edges
        fs << " " << triangle.edges[0] << " " << triangle.edges[1] << " " << triangle.edges[2];
        fs << "\n";
    }
    fs.close();
    return true;
}

std::optional<mesh::TriMesh> Trivis::LoadMesh(
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
        mesh::TriMesh out_mesh;
        while (!fs.eof()) {
            fs >> token;
            if (token == "[NODES]") {
                int n_vertices;
                fs >> n_vertices;
                out_mesh.nodes.resize(n_vertices);
                for (int i = 0; i < n_vertices; ++i) {
                    int node_id;
                    fs >> node_id;
                    auto &node = out_mesh.nodes[node_id];
                    // point
                    fs >> node.point.x >> node.point.y;
                    // edges
                    int size;
                    fs >> size;
                    node.edges.resize(size);
                    for (int &e: node.edges) {
                        fs >> e;
                    }
                    // triangles
                    fs >> size;
                    node.triangles.resize(size);
                    for (int &t: node.triangles) {
                        fs >> t;
                    }
                    // split partner
                    int next_weakly_simple_node;
                    fs >> next_weakly_simple_node;
                    if (next_weakly_simple_node != -1) {
                        node.next_weakly_intersect_node = next_weakly_simple_node;
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
                    // nodes
                    fs >> edge.nodes[0] >> edge.nodes[1];
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
                    // nodes
                    fs >> triangle.nodes[0] >> triangle.nodes[1] >> triangle.nodes[2];
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


///=================///
/// PRIVATE MEMBERS ///
///=================///

/// ##### EXPAND EDGE: INTERSECTION OF RAY AND OBSTACLE ##### ///

Trivis::RayShootingResult Trivis::ExpandEdgeObstacleIntersection(
    int level,
    const geom::FPoint &q,
    const geom::FPoint &distant_t,
    int node_l_id,
    int node_r_id,
    int curr_edge_id,
    int curr_edge_tri_id,
    ExpansionStats *stats
) const {

    if (stats) {
        stats->num_expansions++;
        stats->max_recursion_depth = std::max(stats->max_recursion_depth, level);
    }

    const auto &edge = _mesh.edges[curr_edge_id];
    const auto &node_l_p = _mesh.point(node_l_id);
    const auto &node_r_p = _mesh.point(node_r_id);
    int tri_id = edge.triangles[curr_edge_tri_id];
    const auto &tri = _mesh.triangles[tri_id];
    int opp_id = edge.opposites[curr_edge_tri_id];
    const auto &opp = _mesh.point(opp_id);
    int tri_edge_id = tri.edges[0] == curr_edge_id ? 0 : (tri.edges[1] == curr_edge_id ? 1 : 2);
    // *** note: tri edges are in counter-clockwise order
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    const auto &edge_r = _mesh.edges[edge_r_id];
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    const auto &edge_l = _mesh.edges[edge_l_id];
    if (IsPointInCone(distant_t, node_l_p, q, opp) && TurnsLeft(q, opp, node_l_p)) {
        // tp is in the left cone.
        if (edge_l.is_boundary()) {
            Trivis::RayShootingResult ret;
            ret.code = RaySegmentIntersection(q, distant_t, node_l_p, opp, ret.p, ret.p2);
            if (ret.code == '0' || ret.code == 'c' || ret.code == 'e') {
                return ret;
            }
            ret.edge_id = edge_l_id;
            if (ret.code == 'V' || ret.code == 's' || ret.code == 'v' || ret.code == 'i') {
                if (ret.p == node_l_p) {
                    ret.node_id = node_l_id;
                    return ret;
                }
                if (ret.p == opp) {
                    ret.node_id = opp_id;
                    return ret;
                }
            }
            return ret;
        }
        // expand left edge
        return ExpandEdgeObstacleIntersection(level + 1, q, distant_t, node_l_id, opp_id, edge_l_id, edge_l.triangles[0] == tri_id ? 1 : 0, stats);
    } else { // tp must be in the right cone.
        if (edge_r.is_boundary()) {
            Trivis::RayShootingResult ret;
            ret.code = RaySegmentIntersection(q, distant_t, opp, node_r_p, ret.p, ret.p2);
            if (ret.code == '0' || ret.code == 'c' || ret.code == 'e') {
                return ret;
            }
            ret.edge_id = edge_r_id;
            if (ret.code == 'V' || ret.code == 's' || ret.code == 'v' || ret.code == 'i') {
                if (ret.p == opp) {
                    ret.node_id = opp_id;
                    return ret;
                }
                if (ret.p == node_r_p) {
                    ret.node_id = node_r_id;
                    return ret;
                }
            }
            return ret;
        }
        // expand right edge
        return ExpandEdgeObstacleIntersection(level + 1, q, distant_t, opp_id, node_r_id, edge_r_id, edge_r.triangles[0] == tri_id ? 1 : 0, stats);
    }
}

/// ##### EXPAND EDGE: TWO-POINT QUERIES ##### ///

bool Trivis::ExpandEdgeVisibilityBetween(
    int level,
    const FPoint &q,
    const FPoint &t,
    int node_l_id,
    int node_r_id,
    int curr_edge_id,
    int curr_edge_tri_id,
    ExpansionStats *stats
) const {

    if (stats) {
        stats->num_expansions++;
        stats->max_recursion_depth = std::max(stats->max_recursion_depth, level);
    }

    const auto &edge = _mesh.edges[curr_edge_id];
    const auto &node_l_p = _mesh.point(node_l_id);
    const auto &node_r_p = _mesh.point(node_r_id);
    int tri_id = edge.triangles[curr_edge_tri_id];
    const auto &tri = _mesh.triangles[tri_id];
    int opp_id = edge.opposites[curr_edge_tri_id];
    const auto &opp = _mesh.point(opp_id);
    int tri_edge_id = tri.edges[0] == curr_edge_id ? 0 : (tri.edges[1] == curr_edge_id ? 1 : 2);
    // *** note: tri edges are in counter-clockwise order
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    const auto &edge_r = _mesh.edges[edge_r_id];
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    const auto &edge_l = _mesh.edges[edge_l_id];
    if (IsPointInCone(t, node_l_p, q, opp) && TurnsLeft(q, opp, node_l_p)) {
        // tp is in the left cone.
        if (!TurnsLeft(node_l_p, opp, t)) {
            // tp is in the current triangle.
            return true;
        } else if (!edge_l.is_boundary()) {
            // expand left edge
            return ExpandEdgeVisibilityBetween(level + 1, q, t, node_l_id, opp_id, edge_l_id, edge_l.triangles[0] == tri_id ? 1 : 0, stats);
        } else if (!edge_r.is_boundary() && IsPointInCone(t, opp, q, node_r_p) && TurnsRight(q, opp, node_r_p)) {
            // if the left is obstacle but the right is not and p1 is in both cones, then expand right edge
            return ExpandEdgeVisibilityBetween(level + 1, q, t, opp_id, node_r_id, edge_r_id, edge_r.triangles[0] == tri_id ? 1 : 0, stats);
        } else {
            return false;
        }
    } else { // tp must be in the right cone.
        if (!TurnsRight(node_r_p, opp, t)) {
            // tp is in the current triangle.
            return true;
        } else if (!edge_r.is_boundary()) {
            // expand right edge
            return ExpandEdgeVisibilityBetween(level + 1, q, t, opp_id, node_r_id, edge_r_id, edge_r.triangles[0] == tri_id ? 1 : 0, stats);
        } else {
            return false;
        }
    }
}

/// ##### EXPAND EDGE: MAP VERTICES ##### ///

void Trivis::ExpandEdgeVisibleVertices(
    int level,
    const FPoint &q,
    int rest_l_id,
    int rest_r_id,
    int curr_edge_id,
    int curr_edge_tri_id,
    double sq_radius,
    std::vector<int> &visible_vertices,
    const std::vector<bool> *tabu_vertices,
    ExpansionStats *stats
) const {

    // current edge (struct)
    const auto &curr_edge = _mesh.edges[curr_edge_id];

    bool too_far_away = sq_radius >= 0.0 && PointSegmentSquaredDistance(q, _mesh.point(curr_edge.nodes[0]), _mesh.point(curr_edge.nodes[1])) > sq_radius;

    if (too_far_away || curr_edge.is_boundary()) {
        // Report the edge and return.

        // Find out which one of the edge endpoints is 'left' and which one is 'right'.
        int v_l_id = curr_edge.nodes[0];
        int v_r_id = curr_edge.nodes[1];
        if (_mesh.nodes[v_l_id].edges.front() == curr_edge_id) {
            std::swap(v_l_id, v_r_id);
        }

        if (v_r_id == rest_r_id && (visible_vertices.empty() || visible_vertices.back() != v_r_id) && (sq_radius < 0.0 || q.SquaredDistanceTo(_mesh.point(v_r_id)) <= sq_radius)) {
            if ((!tabu_vertices || !tabu_vertices->operator[](v_r_id)) && (visible_vertices.empty() || visible_vertices.front() != v_r_id)) {
                visible_vertices.push_back(v_r_id);
            }
        }
        if (v_l_id == rest_l_id && (visible_vertices.empty() || visible_vertices.back() != v_l_id) && (sq_radius < 0.0 || q.SquaredDistanceTo(_mesh.point(v_l_id)) <= sq_radius)) {
            if ((!tabu_vertices || !tabu_vertices->operator[](v_l_id)) && (visible_vertices.empty() || visible_vertices.front() != v_l_id)) {
                visible_vertices.push_back(v_l_id);
            }
        }

        return;
    }

    // Record proper expansion (the expansion is not proper if the current edge is an obstacle or too far away).
    if (stats) {
        stats->num_expansions++;
        stats->max_recursion_depth = std::max(stats->max_recursion_depth, level);
    }

    // Expand neighboring edges.

    // current left restriction point
    const auto &rest_l_p = _mesh.point(rest_l_id);
    // current right restriction point
    const auto &rest_r_p = _mesh.point(rest_r_id);

    // current triangle index
    int tri_id = curr_edge.triangles[curr_edge_tri_id];
    // current triangle (struct)
    const auto &tri = _mesh.triangles[tri_id];
    // node (its index) opposite to the current edge in current triangle
    int opp_id = curr_edge.opposites[curr_edge_tri_id];
    // point corresponding to the opposite node
    const auto &opp_p = _mesh.point(opp_id);
    // index (0, 1, or 2) of the current edge for member 'edges' of the current triangle
    int tri_edge_id = tri.edges[0] == curr_edge_id ? 0 : (tri.edges[1] == curr_edge_id ? 1 : 2);

    // index of the edge on right of the current edge
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    // edge (struct) on right of the current edge
    const auto &edge_r = _mesh.edges[edge_r_id];
    // index of the other node (not opp) that makes edge_r
    int edge_r_not_opp_id = edge_r.nodes[edge_r.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_r_is_obstacle = edge_r.is_boundary();

    // index of the edge on left of the current edge
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    // edge (struct) on left of the current edge
    const auto &edge_l = _mesh.edges[edge_l_id];
    // index of the other node (not opp) that makes edge_l
    int edge_l_not_opp_id = edge_l.nodes[edge_l.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_l_is_obstacle = edge_l.is_boundary();

    // true iff the opp_p satisfies restriction given by the left point
    Orientation o_q_rest_l_opp = Orient(q, rest_l_p, opp_p);
    // true iff the opp_p satisfies restriction given by the right point
    Orientation o_q_rest_r_opp = Orient(q, rest_r_p, opp_p);

    // Now decide if the other edges (left=l and right=r of the current edge) of the current triangle will be expanded.

    if (o_q_rest_r_opp != Orientation::kRightTurn) { // the right edge is at least partially visible
        int new_rest_l_id = rest_l_id;
        if (o_q_rest_l_opp != Orientation::kLeftTurn) {
            // update left restriction point if opp is more or the same way restricting
            new_rest_l_id = opp_id;
        }
        int new_edge_tri_id = edge_r.triangles[0] == tri_id ? 1 : 0;

        // *** note: right restriction node will always stay the same

        ExpandEdgeVisibleVertices(level + 1, q, new_rest_l_id, rest_r_id, edge_r_id, new_edge_tri_id, sq_radius, visible_vertices, tabu_vertices, stats);
    }

    if (o_q_rest_l_opp != Orientation::kLeftTurn) { // the left edge is at least partially visible
        int new_rest_r_id = rest_r_id;
        if (o_q_rest_r_opp != Orientation::kRightTurn) {
            // update right restriction point if opp is more or the same way restricting
            new_rest_r_id = opp_id;
        }
        int new_edge_tri_id = edge_l.triangles[0] == tri_id ? 1 : 0;

        // *** note: left restriction node will always stay the same

        ExpandEdgeVisibleVertices(level + 1, q, rest_l_id, new_rest_r_id, edge_l_id, new_edge_tri_id, sq_radius, visible_vertices, tabu_vertices, stats);
    }
}

/// ##### EXPAND EDGE: INPUT POINTS ##### ///

void Trivis::ExpandEdgeVisiblePoints(
    const geom::FPoints &points,
    const std::vector<std::vector<int>> &triangle_points,
    int level,
    const geom::FPoint &q,
    int rest_l_id,
    int rest_r_id,
    int curr_edge_id,
    int curr_edge_tri_id,
    double sq_radius,
    std::vector<bool> &point_visited,
    std::vector<int> &visible_points,
    ExpansionStats *stats
) const {
    // current edge (struct)
    const auto &curr_edge = _mesh.edges[curr_edge_id];

    bool too_far_away = sq_radius >= 0.0 && PointSegmentSquaredDistance(q, _mesh.point(curr_edge.nodes[0]), _mesh.point(curr_edge.nodes[1])) > sq_radius;

    if (too_far_away || curr_edge.is_boundary()) {
        return;
    }

    // Record proper expansion (the expansion is not proper if the current edge is an obstacle or too far away).
    if (stats) {
        stats->num_expansions++;
        stats->max_recursion_depth = std::max(stats->max_recursion_depth, level);
    }

    // Expand neighboring edges.

    // current left restriction point
    const auto &rest_l_p = _mesh.point(rest_l_id);
    // current right restriction point
    const auto &rest_r_p = _mesh.point(rest_r_id);

    // current triangle index
    int tri_id = curr_edge.triangles[curr_edge_tri_id];
    // current triangle (struct)
    const auto &tri = _mesh.triangles[tri_id];
    // node (its index) opposite to the current edge in current triangle
    int opp_id = curr_edge.opposites[curr_edge_tri_id];
    // point corresponding to the opposite node
    const auto &opp_p = _mesh.point(opp_id);
    // index (0, 1, or 2) of the current edge for member 'edges' of the current triangle
    int tri_edge_id = tri.edges[0] == curr_edge_id ? 0 : (tri.edges[1] == curr_edge_id ? 1 : 2);

    // index of the edge on right of the current edge
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    // edge (struct) on right of the current edge
    const auto &edge_r = _mesh.edges[edge_r_id];
    // index of the other node (not opp) that makes edge_r
    int edge_r_not_opp_id = edge_r.nodes[edge_r.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_r_is_obstacle = edge_r.is_boundary();

    // index of the edge on left of the current edge
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    // edge (struct) on left of the current edge
    const auto &edge_l = _mesh.edges[edge_l_id];
    // index of the other node (not opp) that makes edge_l
    int edge_l_not_opp_id = edge_l.nodes[edge_l.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_l_is_obstacle = edge_l.is_boundary();

    // true iff the opp_p satisfies restriction given by the left point
    Orientation o_q_rest_l_opp = Orient(q, rest_l_p, opp_p);
    // true iff the opp_p satisfies restriction given by the right point
    Orientation o_q_rest_r_opp = Orient(q, rest_r_p, opp_p);

    // Report visible points from the current triangle.

    for (int point_id: triangle_points[tri_id]) {
        if (point_visited[point_id]) {
            // Ignore already visited points.
            continue; // Go to the next point
        }
        const auto &p = points[point_id];
        if (sq_radius >= 0.0 && q.SquaredDistanceTo(p) > sq_radius) {
            // The point is too far away: do NOT add, mark as visited.
            point_visited[point_id] = true;
            continue; // Go to the next point
        }
        if (Orient(q, rest_r_p, p) == Orientation::kRightTurn) {
            // The point is to the right of the right restriction ray: do NOT add, do NOT mark as visited. It may be visible by another view.
            continue; // Go to the next point
        }
        if (Orient(q, rest_l_p, p) == Orientation::kLeftTurn) {
            // The point is to the left of the left restriction ray: do NOT add, do NOT mark as visited. It may be visible by another view.
            continue; // Go to the next point
        }
        // The point is in between the restriction rays: add it and mark as visited.
        point_visited[point_id] = true;
        visible_points.push_back(point_id);
        // Go to the next point.
    }

    // Now decide if the other edges (left=l and right=r of the current edge) of the current triangle will be expanded.

    if (o_q_rest_r_opp != Orientation::kRightTurn) { // the right edge is at least partially visible
        int new_rest_l_id = rest_l_id;
        if (o_q_rest_l_opp != Orientation::kLeftTurn) {
            // update left restriction point if opp is more or the same way restricting
            new_rest_l_id = opp_id;
        }
        int new_edge_tri_id = edge_r.triangles[0] == tri_id ? 1 : 0;

        // *** note: right restriction node will always stay the same

        ExpandEdgeVisiblePoints(points, triangle_points, level + 1, q, new_rest_l_id, rest_r_id, edge_r_id, new_edge_tri_id, sq_radius, point_visited, visible_points, stats);
    }

    if (o_q_rest_l_opp != Orientation::kLeftTurn) { // the left edge is at least partially visible
        int new_rest_r_id = rest_r_id;
        if (o_q_rest_r_opp != Orientation::kRightTurn) {
            // update right restriction point if opp is more or the same way restricting
            new_rest_r_id = opp_id;
        }
        int new_edge_tri_id = edge_l.triangles[0] == tri_id ? 1 : 0;

        // *** note: left restriction node will always stay the same

        ExpandEdgeVisiblePoints(points, triangle_points, level + 1, q, rest_l_id, new_rest_r_id, edge_l_id, new_edge_tri_id, sq_radius, point_visited, visible_points, stats);
    }
}

/// ##### EXPAND EDGE: VISIBILITY REGIONS ##### ///

void Trivis::ExpandEdgeVisibilityRegion(
    int level,
    const FPoint &q,
    int rest_l_id,
    int rest_r_id,
    int curr_edge_id,
    int curr_edge_tri_id,
    double sq_radius,
    AbstractVisibilityRegion &visibility_region,
    ExpansionStats *stats
) const {

    // current edge (struct)
    const auto &curr_edge = _mesh.edges[curr_edge_id];

    bool too_far_away = sq_radius >= 0.0 && PointSegmentSquaredDistance(q, _mesh.point(curr_edge.nodes[0]), _mesh.point(curr_edge.nodes[1])) >= sq_radius;

    if (too_far_away || curr_edge.is_boundary()) {
        // Report the edge and return.

        AbstractVisibilityRegionSegment seg;
        seg.id = curr_edge_id;
        if (too_far_away) {
            seg.id = -10 - seg.id;
        }

        // Find out which one of the edge endpoints is 'left' and which one is 'right'.
        int v_l_id = curr_edge.nodes[0];
        int v_r_id = curr_edge.nodes[1];
        if (_mesh.nodes[v_l_id].edges.front() == curr_edge_id) {
            std::swap(v_l_id, v_r_id);
        }

        // Report the first (=right) point.
        if (v_r_id == rest_r_id) {
            seg.v1.id = v_r_id;
        } else {
            // The first point is an intersection of ray (q, right restriction point) and the edge.
            seg.v1.is_intersection = true;
            seg.v1.id_b = rest_r_id;
            seg.v1.id_c = v_l_id;
            seg.v1.id_d = v_r_id;
        }

        // Report the second (=left) point.
        if (v_l_id == rest_l_id) {
            seg.v2.id = v_l_id;
        } else {
            // The second point is an intersection of ray (q, left restriction point) and the edge.
            seg.v2.is_intersection = true;
            seg.v2.id_b = rest_l_id;
            seg.v2.id_c = v_l_id;
            seg.v2.id_d = v_r_id;
        }

        visibility_region.segments.push_back(seg);
        return;
    }

    // Record proper expansion (the expansion is not proper if the current edge is an obstacle or too far away).
    if (stats) {
        stats->num_expansions++;
        stats->max_recursion_depth = std::max(stats->max_recursion_depth, level);
    }

    // Expand neighboring edges.

    // current left restriction point
    const auto &rest_l_p = _mesh.point(rest_l_id);
    // current right restriction point
    const auto &rest_r_p = _mesh.point(rest_r_id);

    // current triangle index
    int tri_id = curr_edge.triangles[curr_edge_tri_id];
    // current triangle (struct)
    const auto &tri = _mesh.triangles[tri_id];
    // node (its index) opposite to the current edge in current triangle
    int opp_id = curr_edge.opposites[curr_edge_tri_id];
    // point corresponding to the opposite node
    const auto &opp_p = _mesh.point(opp_id);
    // index (0, 1, or 2) of the current edge for member 'edges' of the current triangle
    int tri_edge_id = tri.edges[0] == curr_edge_id ? 0 : (tri.edges[1] == curr_edge_id ? 1 : 2);

    // index of the edge on right of the current edge
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    // edge (struct) on right of the current edge
    const auto &edge_r = _mesh.edges[edge_r_id];
    // index of the other node (not opp) that makes edge_r
    int edge_r_not_opp_id = edge_r.nodes[edge_r.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_r_is_obstacle = edge_r.is_boundary();

    // index of the edge on left of the current edge
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    // edge (struct) on left of the current edge
    const auto &edge_l = _mesh.edges[edge_l_id];
    // index of the other node (not opp) that makes edge_l
    int edge_l_not_opp_id = edge_l.nodes[edge_l.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_l_is_obstacle = edge_l.is_boundary();

    // true iff the opp_p satisfies restriction given by the left point
    Orientation o_q_rest_l_opp = Orient(q, rest_l_p, opp_p);
    // true iff the opp_p satisfies restriction given by the right point
    Orientation o_q_rest_r_opp = Orient(q, rest_r_p, opp_p);

    // Now decide if the other edges (left=l and right=r of the current edge) of the current triangle will be expanded.

    if (o_q_rest_r_opp != Orientation::kRightTurn) { // the right edge is at least partially visible
        int new_rest_l_id = rest_l_id;
        if (o_q_rest_l_opp != Orientation::kLeftTurn) {
            // update left restriction point if opp is more or the same way restricting
            new_rest_l_id = opp_id;
        }
        int new_edge_tri_id = edge_r.triangles[0] == tri_id ? 1 : 0;

        // *** note: right restriction node will always stay the same

        ExpandEdgeVisibilityRegion(level + 1, q, new_rest_l_id, rest_r_id, edge_r_id, new_edge_tri_id, sq_radius, visibility_region, stats);
    }

    if (o_q_rest_l_opp != Orientation::kLeftTurn) { // the left edge is at least partially visible
        int new_rest_r_id = rest_r_id;
        if (o_q_rest_r_opp != Orientation::kRightTurn) {
            // update right restriction point if opp is more or the same way restricting
            new_rest_r_id = opp_id;
        }
        int new_edge_tri_id = edge_l.triangles[0] == tri_id ? 1 : 0;

        // *** note: left restriction node will always stay the same

        ExpandEdgeVisibilityRegion(level + 1, q, rest_l_id, new_rest_r_id, edge_l_id, new_edge_tri_id, sq_radius, visibility_region, stats);
    }
}

void Trivis::ExpandEdgeVisibilityRegionIterative(
    const geom::FPoint &q,
    const EdgeExpansionInfo &expansion,
    double sq_radius,
    AbstractVisibilityRegion &visibility_region,
    bool &has_expansion1,
    EdgeExpansionInfo &expansion1,
    bool &has_expansion2,
    EdgeExpansionInfo &expansion2,
    ExpansionStats *stats
) const {

    // current edge (struct)
    const auto &curr_edge = _mesh.edges[expansion.curr_edge_id];

    bool too_far_away = sq_radius >= 0.0 && PointSegmentSquaredDistance(q, _mesh.point(curr_edge.nodes[0]), _mesh.point(curr_edge.nodes[1])) >= sq_radius;

    if (too_far_away || curr_edge.is_boundary()) {
        // Report the edge and return.

        AbstractVisibilityRegionSegment seg;
        seg.id = expansion.curr_edge_id;
        if (too_far_away) {
            seg.id = -10 - seg.id;
        }

        // Find out which one of the edge endpoints is 'left' and which one is 'right'.
        int v_l_id = curr_edge.nodes[0];
        int v_r_id = curr_edge.nodes[1];
        if (_mesh.nodes[v_l_id].edges.front() == expansion.curr_edge_id) {
            std::swap(v_l_id, v_r_id);
        }

        // Report the first (=right) point.
        if (v_r_id == expansion.rest_r_id) {
            seg.v1.id = v_r_id;
        } else {
            // The first point is an intersection of ray (q, right restriction point) and the edge.
            seg.v1.is_intersection = true;
            seg.v1.id_b = expansion.rest_r_id;
            seg.v1.id_c = v_l_id;
            seg.v1.id_d = v_r_id;
        }

        // Report the second (=left) point.
        if (v_l_id == expansion.rest_l_id) {
            seg.v2.id = v_l_id;
        } else {
            // The second point is an intersection of ray (q, left restriction point) and the edge.
            seg.v2.is_intersection = true;
            seg.v2.id_b = expansion.rest_l_id;
            seg.v2.id_c = v_l_id;
            seg.v2.id_d = v_r_id;
        }

        visibility_region.segments.push_back(seg);
        has_expansion1 = false;
        has_expansion2 = false;
        return;
    }

    // Record proper expansion (the expansion is not proper if the current edge is an obstacle or too far away).
    if (stats) {
        stats->num_expansions++;
    }

    // Expand neighboring edges.

    // current left restriction point
    const auto &rest_l_p = _mesh.point(expansion.rest_l_id);
    // current right restriction point
    const auto &rest_r_p = _mesh.point(expansion.rest_r_id);

    // current triangle index
    int tri_id = curr_edge.triangles[expansion.curr_edge_tri_id];
    // current triangle (struct)
    const auto &tri = _mesh.triangles[tri_id];
    // node (its index) opposite to the current edge in current triangle
    int opp_id = curr_edge.opposites[expansion.curr_edge_tri_id];
    // point corresponding to the opposite node
    const auto &opp_p = _mesh.point(opp_id);
    // index (0, 1, or 2) of the current edge for member 'edges' of the current triangle
    int tri_edge_id = tri.edges[0] == expansion.curr_edge_id ? 0 : (tri.edges[1] == expansion.curr_edge_id ? 1 : 2);

    // index of the edge on right of the current edge
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    // edge (struct) on right of the current edge
    const auto &edge_r = _mesh.edges[edge_r_id];
    // index of the other node (not opp) that makes edge_r
    int edge_r_not_opp_id = edge_r.nodes[edge_r.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_r_is_obstacle = edge_r.is_boundary();

    // index of the edge on left of the current edge
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    // edge (struct) on left of the current edge
    const auto &edge_l = _mesh.edges[edge_l_id];
    // index of the other node (not opp) that makes edge_l
    int edge_l_not_opp_id = edge_l.nodes[edge_l.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_l_is_obstacle = edge_l.is_boundary();

    // true iff the opp_p satisfies restriction given by the left point
    Orientation o_q_rest_l_opp = Orient(q, rest_l_p, opp_p);
    // true iff the opp_p satisfies restriction given by the right point
    Orientation o_q_rest_r_opp = Orient(q, rest_r_p, opp_p);

    // Now decide if the other edges (left=l and right=r of the current edge) of the current triangle will be expanded.

    if (o_q_rest_r_opp != Orientation::kRightTurn) { // the right edge is at least partially visible
        int new_rest_l_id = expansion.rest_l_id;
        if (o_q_rest_l_opp != Orientation::kLeftTurn) {
            // update left restriction point if opp is more or the same way restricting
            new_rest_l_id = opp_id;
        }
        int new_edge_tri_id = edge_r.triangles[0] == tri_id ? 1 : 0;

        // *** note: right restriction node will always stay the same

        has_expansion1 = true;
        expansion1.level = expansion.level + 1;
        expansion1.rest_l_id = new_rest_l_id;
        expansion1.rest_r_id = expansion.rest_r_id;
        expansion1.curr_edge_id = edge_r_id;
        expansion1.curr_edge_tri_id = new_edge_tri_id;

    } else {
        has_expansion1 = false;
    }

    if (o_q_rest_l_opp != Orientation::kLeftTurn) { // the left edge is at least partially visible
        int new_rest_r_id = expansion.rest_r_id;
        if (o_q_rest_r_opp != Orientation::kRightTurn) {
            // update right restriction point if opp is more or the same way restricting
            new_rest_r_id = opp_id;
        }
        int new_edge_tri_id = edge_l.triangles[0] == tri_id ? 1 : 0;

        // *** note: left restriction node will always stay the same

        has_expansion2 = true;
        expansion2.level = expansion.level + 1;
        expansion2.rest_l_id = expansion.rest_l_id;
        expansion2.rest_r_id = new_rest_r_id;
        expansion2.curr_edge_id = edge_l_id;
        expansion2.curr_edge_tri_id = new_edge_tri_id;

    } else {
        has_expansion2 = false;
    }
}

void Trivis::ExpandEdgeVisibilityRegionWithHistory(
    int level,
    const FPoint &q,
    int rest_l_id,
    int rest_r_id,
    int curr_edge_id,
    int curr_edge_tri_id,
    double sq_radius,
    AbstractVisibilityRegion &visibility_region,
    std::vector<ExpansionHistoryStep> &history,
    ExpansionStats *stats
) const {

    // current edge (struct)
    const auto &curr_edge = _mesh.edges[curr_edge_id];

    ExpansionHistoryStep history_step;
    history_step.edge_id = curr_edge_id;
    history_step.rest_l_id = rest_l_id;
    history_step.rest_r_id = rest_r_id;

    bool too_far_away = sq_radius >= 0.0 && PointSegmentSquaredDistance(q, _mesh.point(curr_edge.nodes[0]), _mesh.point(curr_edge.nodes[1])) >= sq_radius;

    if (too_far_away || curr_edge.is_boundary()) {
        // Report the edge and return.

        AbstractVisibilityRegionSegment seg;
        seg.id = curr_edge_id;
        if (too_far_away) {
            seg.id = -10 - seg.id;
        }

        // Find out which one of the edge endpoints is 'left' and which one is 'right'.
        int v_l_id = curr_edge.nodes[0];
        int v_r_id = curr_edge.nodes[1];
        if (_mesh.nodes[v_l_id].edges.front() == curr_edge_id) {
            std::swap(v_l_id, v_r_id);
        }

        // Report the first (=right) point.
        if (v_r_id == rest_r_id) {
            seg.v1.id = v_r_id;
        } else {
            // The first point is an intersection of ray (q, right restriction point) and the edge.
            seg.v1.is_intersection = true;
            seg.v1.id_b = rest_r_id;
            seg.v1.id_c = v_l_id;
            seg.v1.id_d = v_r_id;
        }

        // Report the second (=left) point.
        if (v_l_id == rest_l_id) {
            seg.v2.id = v_l_id;
        } else {
            // The second point is an intersection of ray (q, left restriction point) and the edge.
            seg.v2.is_intersection = true;
            seg.v2.id_b = rest_l_id;
            seg.v2.id_c = v_l_id;
            seg.v2.id_d = v_r_id;
        }

        history_step.output_segment_id = static_cast<int>(visibility_region.segments.size());
        history.push_back(history_step);
        visibility_region.segments.push_back(seg);
        return;
    }

    history_step.tri_id = curr_edge.triangles[curr_edge_tri_id];
    history.push_back(history_step);

    // Record proper expansion (the expansion is not proper if the current edge is an obstacle or too far away).
    if (stats) {
        stats->num_expansions++;
        stats->max_recursion_depth = std::max(stats->max_recursion_depth, level);
    }

    // Expand neighboring edges.

    // current left restriction point
    const auto &rest_l_p = _mesh.point(rest_l_id);
    // current right restriction point
    const auto &rest_r_p = _mesh.point(rest_r_id);

    // current triangle index
    int tri_id = curr_edge.triangles[curr_edge_tri_id];
    // current triangle (struct)
    const auto &tri = _mesh.triangles[tri_id];
    // node (its index) opposite to the current edge in current triangle
    int opp_id = curr_edge.opposites[curr_edge_tri_id];
    // point corresponding to the opposite node
    const auto &opp_p = _mesh.point(opp_id);
    // index (0, 1, or 2) of the current edge for member 'edges' of the current triangle
    int tri_edge_id = tri.edges[0] == curr_edge_id ? 0 : (tri.edges[1] == curr_edge_id ? 1 : 2);

    // index of the edge on right of the current edge
    int edge_r_id = tri.edges[(tri_edge_id + 1) % 3];
    // edge (struct) on right of the current edge
    const auto &edge_r = _mesh.edges[edge_r_id];
    // index of the other node (not opp) that makes edge_r
    int edge_r_not_opp_id = edge_r.nodes[edge_r.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_r_is_obstacle = edge_r.is_boundary();

    // index of the edge on left of the current edge
    int edge_l_id = tri.edges[(tri_edge_id + 2) % 3];
    // edge (struct) on left of the current edge
    const auto &edge_l = _mesh.edges[edge_l_id];
    // index of the other node (not opp) that makes edge_l
    int edge_l_not_opp_id = edge_l.nodes[edge_l.nodes[0] == opp_id ? 1 : 0];
    // true if the edge is obstacle
    bool edge_l_is_obstacle = edge_l.is_boundary();

    // true iff the opp_p satisfies restriction given by the left point
    Orientation o_q_rest_l_opp = Orient(q, rest_l_p, opp_p);
    // true iff the opp_p satisfies restriction given by the right point
    Orientation o_q_rest_r_opp = Orient(q, rest_r_p, opp_p);

    // Now decide if the other edges (left=l and right=r of the current edge) of the current triangle will be expanded.

    if (o_q_rest_r_opp != Orientation::kRightTurn) { // the right edge is at least partially visible
        int new_rest_l_id = rest_l_id;
        if (o_q_rest_l_opp != Orientation::kLeftTurn) {
            // update left restriction point if opp is more or the same way restricting
            new_rest_l_id = opp_id;
        }
        int new_edge_tri_id = edge_r.triangles[0] == tri_id ? 1 : 0;

        // *** note: right restriction node will always stay the same

        ExpandEdgeVisibilityRegionWithHistory(level + 1, q, new_rest_l_id, rest_r_id, edge_r_id, new_edge_tri_id, sq_radius, visibility_region, history, stats);
    }

    if (o_q_rest_l_opp != Orientation::kLeftTurn) { // the left edge is at least partially visible
        int new_rest_r_id = rest_r_id;
        if (o_q_rest_r_opp != Orientation::kRightTurn) {
            // update right restriction point if opp is more or the same way restricting
            new_rest_r_id = opp_id;
        }
        int new_edge_tri_id = edge_l.triangles[0] == tri_id ? 1 : 0;

        // *** note: left restriction node will always stay the same

        ExpandEdgeVisibilityRegionWithHistory(level + 1, q, rest_l_id, new_rest_r_id, edge_l_id, new_edge_tri_id, sq_radius, visibility_region, history, stats);
    }
}

/// ##### POSTPROCESSING HELPER METHODS ##### ///

void Trivis::AppendNotCollinearWithIntersection(
    const FPoint &q,
    const AbstractVisibilityRegionVertex &v,
    int edge_flag,
    bool is_last,
    bool line_line_mode,
    RadialVisibilityRegion &visibility_region
) const {
    int vertex_flag;
    FPoint p, b2, aux;
    if (v.is_intersection) {
        // Compute the intersection.
        const auto &b_temp = _mesh.point(v.id_b);
        const auto &b = (q.SquaredDistanceTo(b_temp) < 1e-12) ? b_temp + (b_temp - q).CopyNormalized() : b_temp;
        const auto &c = _mesh.point(v.id_c);
        const auto &d = _mesh.point(v.id_d);
        bool intersection_not_found = false;
        if (line_line_mode) {
            // LineLineIntersectionNotCollinear is faster (but less safe) than RaySegmentIntersection
            intersection_not_found = !LineLineIntersectionNotCollinear(q, b, c, d, p);
        }
        if (!line_line_mode || intersection_not_found) {
            char code = RaySegmentIntersection(q, b, c, d, p, aux);
            intersection_not_found = !(code == '1' || code == 'i' || code == 'v' || code == 's' || code == 'V');
            if (intersection_not_found) {
                std::cerr.precision(std::numeric_limits<double>::max_digits10);
                if (code == 'e') {
                    std::cerr << "[TriVis][AppendNotCollinearWithIntersection] RaySeg intersection NOT found although it should REALLY exist (Code: " << code << "): ";
                } else {
                    std::cerr << "[TriVis][AppendNotCollinearWithIntersection] RaySeg intersection NOT found although it should PROBABLY exist (Code: " << code << "): ";
                }
                std::cerr << std::fixed << q << ", " << b << ", " << c << ", " << d << ".\n";
                return;
            }
        }
        vertex_flag = -1;
    } else {
        p = _mesh.point(v.id);
        vertex_flag = v.id;
    }

    // If the points are exactly equal, ignore the new point.
    if (!visibility_region.vertices.empty()) {
        if (p == visibility_region.vertices.back().point) {
            return;
        }
        if (is_last && p == visibility_region.vertices.front().point) {
            visibility_region.vertices.front().edge_flag = edge_flag;
            return;
        }
    }
    // Not collinear append.
    if (visibility_region.vertices.size() > 1) {
        const auto &p_prev_prev = (visibility_region.vertices.crbegin() + 1)->point;
        const auto &p_prev = visibility_region.vertices.back().point;
        if (edge_flag == visibility_region.vertices.back().edge_flag && Collinear(p_prev_prev, p_prev, p) && Extends(p_prev_prev, p_prev, p)) {
            visibility_region.vertices.back().vertex_flag = vertex_flag;
            visibility_region.vertices.back().point = p;
            return;
        }
    }
    visibility_region.vertices.push_back({vertex_flag, edge_flag, p});
}
