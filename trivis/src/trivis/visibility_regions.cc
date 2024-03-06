/**
 * File:   visibility_region.cc
 *
 * Date:   13.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/visibility_regions.h"

#include <iostream>
#include <algorithm>
#include <cassert>

#include "trivis/geom/robust_geometry.h"
#include "trivis/geom/intersections.h"

using namespace trivis;
using namespace trivis::geom;

bool trivis::IsValid(const RadialVisibilityRegion &region) {
    return !region.vertices.empty() || (region.radius && region.radius >= 0.0);
}

double trivis::Area(
    const RadialVisibilityRegion &region
) {
    assert(IsValid(region));
    double radius_sq = region.radius.value_or(0.0) * region.radius.value_or(0.0);
    int n = static_cast<int>(region.vertices.size());
    if (n <= 1) {
        return M_PI * radius_sq;
    }
    double area = 0;
    for (int i_prev = n - 1, i = 0; i < n; i_prev = i++) {
        const auto &v0 = region.vertices[i_prev];
        const auto &v1 = region.vertices[i];
        area += static_cast<double>(v1.point.x + v0.point.x) * static_cast<double>(v1.point.y - v0.point.y) / 2.0;
        if (v1.edge_flag == -2) {
            double a0 = std::atan2(v0.point.y - region.seed.y, v0.point.x - region.seed.x);
            double a1 = std::atan2(v1.point.y - region.seed.y, v1.point.x - region.seed.x);
            double a_diff = a1 - a0;
            a_diff = std::atan2(std::sin(a_diff), std::cos(a_diff));
            double circular_segment_area = radius_sq * (a_diff - std::sin(a_diff)) / 2.0;
            if (circular_segment_area < 0.0) {
                area += M_PI * radius_sq + circular_segment_area;
            } else {
                area += circular_segment_area;
            }
        }
    }
    return area;
}

void trivis::RemoveAntennas(RadialVisibilityRegion &region) {
    assert(IsValid(region));
    int n_minus_1 = static_cast<int>(region.vertices.size()) - 1;
    std::vector<int> antenna_peaks;
    for (int i_prev = n_minus_1, i = 0; i < region.vertices.size(); i_prev = i++) {
        int i_next = i == n_minus_1 ? 0 : i + 1;
        const auto &v_prev = region.vertices[i_prev];
        const auto &v = region.vertices[i];
        const auto &v_next = region.vertices[i_next];
        if (Collinear(v_prev.point, v.point, v_next.point) && !Extends(v_prev.point, v.point, v_next.point)) {
            antenna_peaks.push_back(i);
        }
    }
    if (antenna_peaks.empty()) {
        return;
    }
    int n_antennas = static_cast<int>(antenna_peaks.size());
    std::vector<int> antenna_1st(n_antennas);
    std::vector<bool> antenna_1st_shifted(n_antennas, false);
    std::vector<int> antenna_end(n_antennas);
    std::vector<bool> antenna_end_shifted(n_antennas, false);
    for (int i = 0; i < n_antennas; ++i) {
        int antenna_peak = antenna_peaks[i];
        int antenna_1st_candid, antenna_end_candid;
        {   // Compute antenna_1st_candid.
            int j_curr = antenna_peak;
            while (true) {
                int j_prev = j_curr == 0 ? n_minus_1 : j_curr - 1;
                int j_prev_prev = j_prev == 0 ? n_minus_1 : j_prev - 1;
                const auto &v_curr = region.vertices[j_curr];
                const auto &v_prev = region.vertices[j_prev];
                const auto &v_prev_prev = region.vertices[j_prev_prev];
                if (!Collinear(v_curr.point, v_prev.point, v_prev_prev.point)) {
                    antenna_1st_candid = j_prev;
                    break;
                }
                j_curr = j_prev;
            }
        }
        {   // Compute antenna_end_candid.
            int j_curr = antenna_peak;
            while (true) {
                int j_next = j_curr == n_minus_1 ? 0 : j_curr + 1;
                int j_next_next = j_next == n_minus_1 ? 0 : j_next + 1;
                const auto &v_curr = region.vertices[j_curr];
                const auto &v_next = region.vertices[j_next];
                const auto &v_next_next = region.vertices[j_next_next];
                if (!Collinear(v_curr.point, v_next.point, v_next_next.point)) {
                    antenna_end_candid = j_next;
                    break;
                }
                j_curr = j_next;
            }
        }
        {   // Get antenna real 1st.
            int j_curr = antenna_1st_candid;
            while (j_curr < antenna_peak) {
                int j_next = j_curr == n_minus_1 ? 0 : j_curr + 1;
                const auto &v_curr = region.vertices[j_curr];
                const auto &v_next = region.vertices[j_next];
                if (!Extends(v_curr.point, v_next.point, region.vertices[antenna_end_candid].point)) {
                    antenna_1st[i] = j_curr;
                    break;
                }
                j_curr = j_next;
            }

        }
        antenna_1st_shifted[i] = (antenna_1st_candid != antenna_1st[i]);
        {   // Get antenna real end.
            int j_curr = antenna_end_candid;
            while (j_curr > antenna_peak) {
                int j_prev = j_curr == 0 ? n_minus_1 : j_curr - 1;
                const auto &v_curr = region.vertices[j_curr];
                const auto &v_prev = region.vertices[j_prev];
                if (!Extends(v_curr.point, v_prev.point, region.vertices[antenna_1st_candid].point)) {
                    antenna_end[i] = j_curr;
                    break;
                }
                j_curr = j_prev;
            }

        }
        antenna_end_shifted[i] = (antenna_end_candid != antenna_end[i]);
    }

    int antenna_cnt = 0;
    bool is_antenna = false;
    std::vector<VisibilityRegionVertex> new_vertices;
    for (int i = 0; i < region.vertices.size(); ++i) {
        const auto &v = region.vertices[i];
        if (is_antenna) {
            if (i == antenna_end[antenna_cnt]) {
                new_vertices.push_back(v);
                if (antenna_1st_shifted[antenna_cnt]) {
                    new_vertices.back().edge_flag = region.vertices[antenna_1st[antenna_cnt] == n_minus_1 ? 0 : antenna_1st[antenna_cnt] + 1].edge_flag;
                }
                ++antenna_cnt;
                is_antenna = false;
                continue;
            }
            if (i == n_minus_1) {
                std::reverse(new_vertices.begin(), new_vertices.end());
                for (int j = 0; j < antenna_end[antenna_cnt]; ++j) {
                    new_vertices.pop_back();
                }
                if (antenna_1st_shifted[antenna_cnt]) {
                    new_vertices.back().edge_flag = region.vertices[antenna_1st[antenna_cnt] == n_minus_1 ? 0 : antenna_1st[antenna_cnt] + 1].edge_flag;
                }
                std::reverse(new_vertices.begin(), new_vertices.end());
                break;
            }
            continue;
        }
        if (antenna_cnt < n_antennas && i == antenna_1st[antenna_cnt]) {
            new_vertices.push_back(v);
            is_antenna = true;
            continue;
        }
        new_vertices.push_back(v);
    }
    region.vertices = std::move(new_vertices);
}

void trivis::IntersectWithCircleCenteredAtSeed(
    RadialVisibilityRegion &region,
    std::optional<double> radius
) {
    assert(IsValid(region));
    assert((region.radius && region.radius >= 0.0) || (radius && radius >= 0.0));
    double radius_final = radius.value_or(region.radius.value_or(0.0));
    double sq_radius = radius_final * radius_final;
    int n_minus_1 = static_cast<int>(region.vertices.size()) - 1;
    std::vector<VisibilityRegionVertex> new_vertices;
    for (int i_prev = n_minus_1, i = 0; i < region.vertices.size(); i_prev = i++) {
        const auto &vi_prev = region.vertices[i_prev];
        const auto &vi = region.vertices[i];
        if (vi.edge_flag < -1) {
            // The whole edge is outside the radius: IGNORE IT.
            continue;
        }
        bool is_inside_pi_prev = vi_prev.point.SquaredDistanceTo(region.seed) <= sq_radius;
        bool is_inside_pi = vi.point.SquaredDistanceTo(region.seed) <= sq_radius;
        if (is_inside_pi_prev) {
            if (is_inside_pi) {
                // The whole edge is inside the radius (no intersection).
                if (vi_prev.edge_flag < -1) {
                    new_vertices.push_back({vi_prev.vertex_flag, -2, vi_prev.point});
                }
                new_vertices.push_back({vi.vertex_flag, vi.edge_flag, vi.point});
            } else {
                // The first endpoint is inside, the second is outside (1 intersection).
                FPoint intersection;
                if (vi.edge_flag < 0) {
                    intersection = vi_prev.point - region.seed;
                    intersection = intersection / intersection.Norm();
                    intersection = region.seed + intersection * radius_final;
                } else {
                    auto intersections = LineCircleIntersections(vi_prev.point, vi.point, region.seed, radius_final, true);
                    if (intersections.empty()) {
                        intersections = LineCircleIntersections(vi_prev.point.CopySwappedXY(), vi.point.CopySwappedXY(), region.seed.CopySwappedXY(), radius_final, true);
                        if (intersections.empty()) {
                            std::cerr << std::fixed << "[TriVis][IntersectWithCircle] LineCircle intersections should not be empty (1)! ";
                            std::cerr << "Line: " << vi_prev.point << ", " << vi.point << ", Circle:" << region.seed << ", " << radius_final << ".\n";
                            continue;
                        }
                        intersections[0].SwapXY();
                    }
                    intersection = intersections[0];
                }
                new_vertices.push_back({-1, vi.edge_flag, intersection});
            }
        } else if (is_inside_pi) {
            // The first endpoint is outside, the second is inside (1 intersection).
            FPoint intersection;
            if (vi.edge_flag < 0) {
                intersection = vi_prev.point - region.seed;
                intersection = intersection / intersection.Norm();
                intersection = region.seed + intersection * radius_final;
            } else {
                auto intersections = LineCircleIntersections(vi_prev.point, vi.point, region.seed, radius_final, true);
                if (intersections.empty()) {
                    intersections = LineCircleIntersections(vi_prev.point.CopySwappedXY(), vi.point.CopySwappedXY(), region.seed.CopySwappedXY(), radius_final, true);
                    if (intersections.empty()) {
                        std::cerr << std::fixed << "[TriVis][IntersectWithCircle] LineCircle intersections should not be empty (2)! ";
                        std::cerr << "Line: " << vi_prev.point << ", " << vi.point << ", Circle:" << region.seed << ", " << radius_final << ".\n";
                        continue;
                    }
                    intersections[0].SwapXY();
                }
                intersection = intersections[0];
            }
            new_vertices.push_back({-1, -2, intersection});
            new_vertices.push_back(vi);
        } else {
            // Both endpoints are outside (2 intersections).
            auto intersections = LineCircleIntersections(vi_prev.point, vi.point, region.seed, radius_final, true);
            if (!intersections.empty()) {
                if (intersections.size() > 1 && vi_prev.point.SquaredDistanceTo(intersections[0]) > vi_prev.point.SquaredDistanceTo(intersections[1])) {
                    std::swap(intersections[0], intersections[1]);
                }
                new_vertices.push_back({-1, -2, intersections[0]});
                if (intersections.size() > 1) {
                    new_vertices.push_back({-1, vi.edge_flag, intersections[1]});
                }
            }
        }
    }
    region.radius = radius_final;
    region.vertices = std::move(new_vertices);
}

void trivis::RemoveShortEdges(
    RadialVisibilityRegion &region,
    double min_edge_length
) {
    assert(IsValid(region));
    std::vector<VisibilityRegionVertex> new_vertices;
    for (int i_prev = static_cast<int>(region.vertices.size()) - 1, i = 0; i < region.vertices.size(); i_prev = i++) {
        const auto &v_prev = region.vertices[i_prev];
        const auto &v = region.vertices[i];
        if (v_prev.point.DistanceTo(v.point) < min_edge_length) {
            if (v.vertex_flag >= 0) {
                if (!new_vertices.empty()) {
                    new_vertices.pop_back();
                }
                new_vertices.push_back({v.vertex_flag, v_prev.edge_flag, v.point});
            }
        } else {
            new_vertices.push_back(v);
        }
    }
    region.vertices = std::move(new_vertices);
}

void trivis::SampleArcEdges(
    RadialVisibilityRegion &region,
    double max_sample_angle
) {
    assert(IsValid(region));
    if (!region.radius) {
        return;
    }
    double radius = *region.radius;
    std::vector<VisibilityRegionVertex> new_vertices;
    int n_v = static_cast<int>(region.vertices.size());
    if (n_v <= 1) {
        double a_diff = 2 * M_PI;
        int angle_n = static_cast<int>(std::ceil(a_diff / max_sample_angle));
        double angle_delta = a_diff / static_cast<double>(angle_n);
        double alpha = 0;
        for (int j = 0; j < angle_n - 1; ++j) {
            alpha += angle_delta;
            FPoint p = MakePoint(region.seed.x + radius * std::cos(alpha), region.seed.y + radius * std::sin(alpha));
            new_vertices.push_back({-1, -4, p});
        }
        region.vertices = std::move(new_vertices);
        return;
    }
    for (int i_prev = n_v - 1, i = 0; i < n_v; i_prev = i++) {
        const auto v_prev = region.vertices[i_prev];
        const auto v = region.vertices[i];
        if (v.edge_flag != -2) {
            new_vertices.push_back(v);
            continue;
        }
        double a0 = std::atan2(v_prev.point.y - region.seed.y, v_prev.point.x - region.seed.x);
        double a1 = std::atan2(v.point.y - region.seed.y, v.point.x - region.seed.x);
        if (a1 < a0) {
            a1 += 2.0 * M_PI;
        }
        double a_diff = a1 - a0;
        int angle_n = static_cast<int>(std::ceil(a_diff / max_sample_angle));
        double angle_delta = a_diff / static_cast<double>(angle_n);
        double alpha = a0;
        for (int j = 0; j < angle_n - 1; ++j) {
            alpha += angle_delta;
            FPoint p = MakePoint(region.seed.x + radius * std::cos(alpha), region.seed.y + radius * std::sin(alpha));
            new_vertices.push_back({-1, -4, p});
        }
        new_vertices.push_back(v);
        new_vertices.back().edge_flag = -4;
    }
    region.vertices = std::move(new_vertices);
}

geom::FPolygon trivis::ToPolygon(
    const RadialVisibilityRegion &region
) {
    assert(IsValid(region));
    geom::FPolygon ret;
    ret.reserve(region.vertices.size());
    for (const auto &v: region.vertices) {
        ret.push_back(v.point);
    }
    return ret;
}