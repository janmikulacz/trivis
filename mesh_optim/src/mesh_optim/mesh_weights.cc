/**
 * File:   mesh_weights.cc
 *
 * Date:   20.10.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "mesh_optim/mesh_weights.h"

#include <thread>
#include <cassert>

#include "trivis/utils/clipper_utils.h"

using namespace trivis;
using namespace trivis::geom;
using namespace mesh_optim;

std::vector<std::vector<double>> mesh_optim::ComputeWeightsEdgeLength(
    const Trivis &vis,
    double eps_dist_diff_collinear
) {
    int n_v = static_cast<int>(vis.mesh().vertices.size());

    auto weights = std::vector<std::vector<double>>(n_v, std::vector<double>(n_v, kMaxDouble));
    auto dists = std::vector<std::vector<double>>(n_v, std::vector<double>(n_v, kMaxDouble));

    // Compute the weights.
    for (int i = 0; i < n_v; ++i) {
        const auto &pi = vis.mesh().point(i);
        std::vector<int> visible_vertices = vis.VisibleVertices(i);
        for (int j: visible_vertices) {
            if (weights[i][j] == kMaxDouble) {
                if (i == j) {
                    weights[i][j] = weights[j][i] = 0.0;
                    dists[i][j] = dists[j][i] = 0.0;
                    continue;
                }
                bool is_edge_obstacle = false;
                for (int e_idx: vis.mesh().vertices[i].edges) {
                    const auto &e = vis.mesh().edges[e_idx];
                    if (e.is_boundary() && ((e.vertices[0] == i && e.vertices[1] == j) || (e.vertices[0] == j && e.vertices[1] == i))) {
                        is_edge_obstacle = true;
                        break;
                    }
                }
                const auto &pj = vis.mesh().point(j);
                double dist = pi.DistanceTo(pj);
                if (is_edge_obstacle) {
                    weights[i][j] = weights[j][i] = 0.0;
                } else {
                    weights[i][j] = weights[j][i] = dist;
                }
                dists[i][j] = dists[j][i] = dist;
            }
        }
    }

    // Brute force check for close-to-collinear edges.
    for (int i = 0; i < n_v; ++i) {
        for (int j = i + 1; j < n_v; ++j) {
            if (dists[i][j] == kMaxDouble) {
                continue;
            }
            for (int k = 0; k < n_v; ++k) {
                if (i == k || j == k || dists[i][k] == kMaxDouble || dists[j][k] == kMaxDouble) {
                    continue;
                }
                if (std::abs(dists[i][k] + dists[j][k] - dists[i][j]) < eps_dist_diff_collinear) {
                    weights[i][j] = weights[j][i] = kMaxDouble;
                }
            }
        }
    }

    return weights;
}

void ComputeEdgeVisibleArea(
    const Trivis &vis,
    int i,
    int j,
    std::optional<bool> is_obstacle_opt,
    double total_area,
    double max_sample_dist,
    std::optional<double> vis_radius_opt,
    std::vector<std::vector<double>> &weights
) {
    bool is_edge_obstacle = false;
    if (is_obstacle_opt) {
        is_edge_obstacle = *is_obstacle_opt;
    } else {
        for (int e_idx: vis.mesh().vertices[i].edges) {
            const auto &e = vis.mesh().edges[e_idx];
            if (e.is_boundary() && ((e.vertices[0] == i && e.vertices[1] == j) || (e.vertices[0] == j && e.vertices[1] == i))) {
                is_edge_obstacle = true;
                break;
            }
        }
    }
    if (is_edge_obstacle) {
        weights[i][j] = weights[j][i] = 0.0;
    } else {
        const auto &pi = vis.mesh().point(i);
        const auto &pj = vis.mesh().point(j);
        auto sample_points = SampleSegment2D(pi, pj, max_sample_dist, false, false);
        if (sample_points.empty()) {
            sample_points = {(pi + pj) / 2.0};
        }
        int n_sample_points = static_cast<int>(sample_points.size());
        Clipper2Lib::Paths64 solution;
        Clipper2Lib::Paths64 vis_polygons64;
        Clipper2Lib::Clipper64 clipper64;
        for (int k = 0; k < n_sample_points; ++k) {
            const auto &sample_p = sample_points[k];
            auto p_locate = vis.LocatePoint(sample_p);
            if (p_locate) {
                auto abs_reg_opt = vis.VisibilityRegion(sample_p, p_locate.value(), vis_radius_opt);
                vis_polygons64.push_back(Clipper2Lib::PathDToPath64(utils::Geom2Clipper(vis.ToRadialVisibilityRegion(abs_reg_opt).ToPolygon())));
                if (vis_polygons64.size() == 10 || k + 1 == n_sample_points) {
                    clipper64.AddClip(solution);
                    clipper64.AddSubject(vis_polygons64);
                    clipper64.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero, solution);
                    clipper64.Clear();
                    vis_polygons64.clear();
                }
            }
        }
        if (vis_radius_opt) {
            geom::Polygons<double> clips;
            clips.push_back(geom::MakeRegularPolygon(pi, *vis_radius_opt, 100)); // fixme: this should be a real parameter
            auto u = pj - pi;
            u = (u / u.Norm()) * (*vis_radius_opt);
            std::swap(u.x, u.y);
            u.x = -u.x;
            clips.push_back(geom::Polygon<double>{pi + u, pj + u, pj - u, pi - u});
            clips.push_back(geom::MakeRegularPolygon(pj, *vis_radius_opt, 100)); // fixme: this should be a real parameter
            clipper64.AddSubject(Clipper2Lib::PathsDToPaths64(utils::Geom2Clipper(clips)));
            Clipper2Lib::Paths64 clip;
            clipper64.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero, clip);
            clipper64.Clear();
            clipper64.AddSubject(solution);
            clipper64.AddClip(clip);
            clipper64.Execute(Clipper2Lib::ClipType::Intersection, Clipper2Lib::FillRule::NonZero, solution);
            clipper64.Clear();
        }
        weights[i][j] = weights[j][i] = Clipper2Lib::Area(solution) / total_area;
    }
}

std::vector<std::vector<double>> mesh_optim::ComputeWeightsEdgeVisibility(
    const Trivis &vis,
    double max_sample_dist,
    double eps_dist_diff_collinear,
    std::optional<double> vis_radius,
    bool parallelize
) {
    int n_v = static_cast<int>(vis.mesh().vertices.size());

    auto weights = std::vector<std::vector<double>>(n_v, std::vector<double>(n_v, kMaxDouble));
    auto dists = std::vector<std::vector<double>>(n_v, std::vector<double>(n_v, kMaxDouble));

    double total_area = Clipper2Lib::Area(Clipper2Lib::PathsDToPaths64(utils::Geom2Clipper(vis.triangles())));

    // Compute the weights.
    std::vector<std::thread> threads;
    for (int i = 0; i < n_v; ++i) {
        const auto &pi = vis.mesh().point(i);
        std::vector<int> visible_vertices = vis.VisibleVertices(i);
        for (int j: visible_vertices) {
            if (dists[i][j] == kMaxDouble) {
                if (i == j) {
                    weights[i][j] = weights[j][i] = 0.0;
                    dists[i][j] = dists[j][i] = 0.0;
                    continue;
                }
                const auto &pj = vis.mesh().point(j);
                double dist = pi.DistanceTo(pj);

                if (parallelize) {
                    threads.emplace_back(ComputeEdgeVisibleArea, std::ref(vis), i, j, std::nullopt, total_area, max_sample_dist, vis_radius, std::ref(weights));
                } else {
                    ComputeEdgeVisibleArea(vis, i, j, std::nullopt, total_area, max_sample_dist, vis_radius, weights);
                }

                dists[i][j] = dists[j][i] = dist;
            }
        }
    }
    if (parallelize) {
        for (auto &thread: threads) {
            thread.join();
        }
    }

    // Brute force check for close-to-collinear edges.
    for (int i = 0; i < n_v; ++i) {
        for (int j = i + 1; j < n_v; ++j) {
            if (dists[i][j] == kMaxDouble) {
                continue;
            }
            for (int k = 0; k < n_v; ++k) {
                if (i == k || j == k || dists[i][k] == kMaxDouble || dists[j][k] == kMaxDouble) {
                    continue;
                }
                if (std::abs(dists[i][k] + dists[j][k] - dists[i][j]) < eps_dist_diff_collinear) {
                    weights[i][j] = weights[j][i] = kMaxDouble;
                }
            }
        }
    }

    return weights;
}

std::vector<std::vector<double>> mesh_optim::ComputeWeightsEdgeVisibilityOptimized(
    const Trivis &vis,
    double max_sample_dist,
    double eps_dist_diff_collinear,
    double long_edge_penalty,
    std::optional<double> long_edge_threshold_opt,
    std::optional<int> shortest_edges_ver_max_opt,
    std::optional<double> vis_radius_opt
) {
    assert(long_edge_threshold_opt || shortest_edges_ver_max_opt);
    assert(!long_edge_threshold_opt || (0.0 <= long_edge_threshold_opt && long_edge_threshold_opt <= 1.0));
    assert(!shortest_edges_ver_max_opt || shortest_edges_ver_max_opt > 0);

    int n_v = static_cast<int>(vis.mesh().vertices.size());

    auto weights = std::vector<std::vector<double>>(n_v, std::vector<double>(n_v, kMaxDouble));
    auto lengths = std::vector<std::vector<double>>(n_v, std::vector<double>(n_v, kMaxDouble));

    double total_area = Clipper2Lib::Area(Clipper2Lib::PathsDToPaths64(utils::Geom2Clipper(vis.triangles())));

    // Precompute edges.
    struct EdgeInfo {
        int i;
        int j;
        double length;
    };
    std::vector<EdgeInfo> free_edges;
    for (int ver_id = 0; ver_id < n_v; ++ver_id) {
        const auto &ver = vis.mesh().vertices[ver_id];
        std::vector<int> obstacle_edges;
        for (int e_id: ver.edges) {
            if (vis.mesh().edges[e_id].is_boundary()) {
                obstacle_edges.push_back(e_id);
            }
        }
        std::vector<int> visible_vertices = vis.VisibleVertices(ver_id);
        std::vector<int> free_visible_vertices;
        std::vector<std::vector<bool>> in_free_edges;
        if (shortest_edges_ver_max_opt) {
            in_free_edges = std::vector<std::vector<bool>>(n_v, std::vector<bool>(n_v, false));
        }
        for (int visible_ver_id: visible_vertices) {
            const auto &visible_ver = vis.mesh().vertices[visible_ver_id];
            if (lengths[ver_id][visible_ver_id] != kMaxDouble) {
                continue;
            }
            double length = ver.point.DistanceTo(visible_ver.point);
            lengths[ver_id][visible_ver_id] = lengths[visible_ver_id][ver_id] = length;
            bool is_obstacle = false;
            for (int e_id: obstacle_edges) {
                const auto &e = vis.mesh().edges[e_id];
                if (e.vertices[e.vertices[0] == ver_id ? 1 : 0] == visible_ver_id) {
                    is_obstacle = true;
                    break;
                }
            }
            if (is_obstacle) {
                // Set obstacle edges weights to zero.
                weights[ver_id][visible_ver_id] = weights[visible_ver_id][ver_id] = 0.0;
            } else if (shortest_edges_ver_max_opt) {
                free_visible_vertices.push_back(visible_ver_id);
            } else {
                free_edges.push_back(EdgeInfo{ver_id, visible_ver_id, length});
            }
        }
        if (shortest_edges_ver_max_opt) {
            const auto &ver_lengths = lengths[ver_id];
            int n_free_visible_vertices = static_cast<int>(free_visible_vertices.size());
            int n_shortest_edges_ver = std::min(*shortest_edges_ver_max_opt, n_free_visible_vertices);
            std::partial_sort(free_visible_vertices.begin(), free_visible_vertices.begin() + n_shortest_edges_ver, free_visible_vertices.end(), [ver_lengths](
                int ver_id_0, int ver_id_1
            ) {
                return ver_lengths[ver_id_0] < ver_lengths[ver_id_1];
            });
            for (int i = 0; i < n_shortest_edges_ver; ++i) {
                int visible_ver_id = free_visible_vertices[i];
                if (!in_free_edges[ver_id][visible_ver_id]) {
                    in_free_edges[ver_id][visible_ver_id] = in_free_edges[visible_ver_id][ver_id] = true;
                    free_edges.push_back(EdgeInfo{ver_id, visible_ver_id, ver_lengths[visible_ver_id]});
                }
            }
            for (int i = n_shortest_edges_ver; i < n_free_visible_vertices; ++i) {
                int visible_ver_id = free_visible_vertices[i];
                if (!in_free_edges[ver_id][visible_ver_id]) {
                    weights[ver_id][visible_ver_id] = weights[visible_ver_id][ver_id] = long_edge_penalty + ver_lengths[visible_ver_id];
                }
            }
        }
    }

    int n_free_edges = static_cast<int>(free_edges.size());
    if (long_edge_threshold_opt) {
        int n_sorted_edges = static_cast<int>(std::ceil(*long_edge_threshold_opt * static_cast<double>(n_free_edges)));
        std::partial_sort(free_edges.begin(), free_edges.begin() + n_sorted_edges, free_edges.end(), [](
            const EdgeInfo &e0, const EdgeInfo &e1
        ) {
            return e0.length < e1.length;
        });
        // Compute the weights of free edges.
        for (int i = 0; i < n_sorted_edges; ++i) {
            const auto &edge = free_edges[i];
            ComputeEdgeVisibleArea(vis, edge.i, edge.j, false, total_area, max_sample_dist, vis_radius_opt, weights);
        }
        // Set penalties for too long edges.
        for (int i = n_sorted_edges; i < n_free_edges; ++i) {
            const auto &edge = free_edges[i];
            weights[edge.i][edge.j] = weights[edge.j][edge.i] = long_edge_penalty + edge.length;
        }
    } else {
        // Compute the weights of free edges.
        for (int i = 0; i < n_free_edges; ++i) {
            const auto &edge = free_edges[i];
            ComputeEdgeVisibleArea(vis, edge.i, edge.j, false, total_area, max_sample_dist, vis_radius_opt, weights);
        }
    }

    // Brute force check for close-to-collinear edges.
    for (int i = 0; i < n_v; ++i) {
        for (int j = i + 1; j < n_v; ++j) {
            if (weights[i][j] == kMaxDouble) {
                continue;
            }
            for (int k = 0; k < n_v; ++k) {
                if (i == k || j == k || lengths[i][k] == kMaxDouble || lengths[j][k] == kMaxDouble) {
                    continue;
                }
                if (std::abs(lengths[i][k] + lengths[j][k] - lengths[i][j]) < eps_dist_diff_collinear) {
                    weights[i][j] = weights[j][i] = kMaxDouble;
                }
            }
        }
    }

    return weights;
}

void mesh_optim::NegateWeights(
    std::vector<std::vector<double>> &weights
) {
    for (auto &weight_row: weights) {
        for (double &c: weight_row) {
            if (c != kMaxDouble) {
                c = -c;
            }
        }
    }
}