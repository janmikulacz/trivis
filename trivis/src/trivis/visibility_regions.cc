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

#include "trivis/geom/robust_geometry.h"

using namespace trivis;
using namespace trivis::geom;

bool trivis::operator==(const AbstractVisibilityRegionVertex &v1, const AbstractVisibilityRegionVertex &v2) {
    bool is_intersection = v1.is_intersection;
    if (v2.is_intersection != is_intersection) {
        // The intersection flag must be the same.
        return false;
    }
    if (!is_intersection) {
        // If not intersection, then id_b, id_c, id_d does not matter.
        return v1.id == v2.id;
    }
    return v1.id == v2.id && v1.id_b == v2.id_b && v1.id_c == v2.id_c && v1.id_d == v2.id_d;
}

bool trivis::operator!=(const AbstractVisibilityRegionVertex &v1, const AbstractVisibilityRegionVertex &v2) {
    return !(v1 == v2);
}

bool trivis::operator==(const AbstractVisibilityRegionSegment &s1, const AbstractVisibilityRegionSegment &s2) {
    return s1.id == s2.id && s1.v1 == s2.v1 && s1.v2 == s2.v2;
}

bool trivis::operator!=(const AbstractVisibilityRegionSegment &s1, const AbstractVisibilityRegionSegment &s2) {
    return !(s1 == s2);
}

bool trivis::operator==(const AbstractVisibilityRegion &r1, const AbstractVisibilityRegion &r2) {
    return r1.seed_id == r2.seed_id && r1.seed == r2.seed && r1.segments == r2.segments;
}

bool trivis::operator!=(const AbstractVisibilityRegion &r1, const AbstractVisibilityRegion &r2) {
    return !(r1 == r2);
}

bool trivis::IsValid(const RadialVisibilityRegion &reg) {
    return !reg.vertices.empty() || reg.radius > 0.0;
}

void trivis::RemoveAntennas(
    RadialVisibilityRegion &res
) {
    int n_minus_1 = static_cast<int>(res.vertices.size()) - 1;
    std::vector<int> antenna_peaks;
    for (int i_prev = n_minus_1, i = 0; i < res.vertices.size(); i_prev = i++) {
        int i_next = i == n_minus_1 ? 0 : i + 1;
        const auto &v_prev = res.vertices[i_prev];
        const auto &v = res.vertices[i];
        const auto &v_next = res.vertices[i_next];
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
                const auto &v_curr = res.vertices[j_curr];
                const auto &v_prev = res.vertices[j_prev];
                const auto &v_prev_prev = res.vertices[j_prev_prev];
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
                const auto &v_curr = res.vertices[j_curr];
                const auto &v_next = res.vertices[j_next];
                const auto &v_next_next = res.vertices[j_next_next];
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
                const auto &v_curr = res.vertices[j_curr];
                const auto &v_next = res.vertices[j_next];
                if (!Extends(v_curr.point, v_next.point, res.vertices[antenna_end_candid].point)) {
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
                const auto &v_curr = res.vertices[j_curr];
                const auto &v_prev = res.vertices[j_prev];
                if (!Extends(v_curr.point, v_prev.point, res.vertices[antenna_1st_candid].point)) {
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
    for (int i = 0; i < res.vertices.size(); ++i) {
        const auto &v = res.vertices[i];
        if (is_antenna) {
            if (i == antenna_end[antenna_cnt]) {
                new_vertices.push_back(v);
                if (antenna_1st_shifted[antenna_cnt]) {
                    new_vertices.back().edge_flag = res.vertices[antenna_1st[antenna_cnt] == n_minus_1 ? 0 : antenna_1st[antenna_cnt] + 1].edge_flag;
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
                    new_vertices.back().edge_flag = res.vertices[antenna_1st[antenna_cnt] == n_minus_1 ? 0 : antenna_1st[antenna_cnt] + 1].edge_flag;
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
    res.vertices = std::move(new_vertices);
}

double trivis::Area(
    const RadialVisibilityRegion &vis_reg
) {
    double radius_sq = vis_reg.radius * vis_reg.radius;
    int n = static_cast<int>(vis_reg.vertices.size());
    if (n <= 1) {
        return M_PI * radius_sq;
    }
    double area = 0;
    for (int i_prev = n - 1, i = 0; i < n; i_prev = i++) {
        const auto &v0 = vis_reg.vertices[i_prev];
        const auto &v1 = vis_reg.vertices[i];
        area += static_cast<double>(v1.point.x + v0.point.x) * static_cast<double>(v1.point.y - v0.point.y) / 2.0;
        if (v1.edge_flag == -2) {
            double a0 = std::atan2(v0.point.y - vis_reg.seed.y, v0.point.x - vis_reg.seed.x);
            double a1 = std::atan2(v1.point.y - vis_reg.seed.y, v1.point.x - vis_reg.seed.x);
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

void trivis::RemoveShortEdges(
    double min_edge_length,
    RadialVisibilityRegion &res
) {
    std::vector<VisibilityRegionVertex> new_vertices;
    for (int i_prev = static_cast<int>(res.vertices.size()) - 1, i = 0; i < res.vertices.size(); i_prev = i++) {
        const auto &v_prev = res.vertices[i_prev];
        const auto &v = res.vertices[i];
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
    res.vertices = std::move(new_vertices);
}

RadialVisibilityRegion trivis::SampleArcEdges(
    const RadialVisibilityRegion &reg,
    double max_sample_beta
) {
    RadialVisibilityRegion ret;
    ret.seed_id = reg.seed_id;
    ret.seed = reg.seed;
    ret.radius = reg.radius;
    int n_v = static_cast<int>(reg.vertices.size());
    if (n_v <= 1) {
        double a_diff = 2 * M_PI;
        int angle_n = static_cast<int>(std::ceil(a_diff / max_sample_beta));
        double angle_delta = a_diff / static_cast<double>(angle_n);
        double alpha = 0;
        for (int j = 0; j < angle_n - 1; ++j) {
            alpha += angle_delta;
            FPoint p = MakePoint(reg.seed.x + reg.radius * std::cos(alpha), reg.seed.y + reg.radius * std::sin(alpha));
            ret.vertices.push_back({-1, -4, p});
        }
        return ret;
    }
    for (int i_prev = n_v - 1, i = 0; i < n_v; i_prev = i++) {
        const auto v_prev = reg.vertices[i_prev];
        const auto v = reg.vertices[i];
        if (v.edge_flag != -2) {
            ret.vertices.push_back(v);
            continue;
        }
        double a0 = std::atan2(v_prev.point.y - reg.seed.y, v_prev.point.x - reg.seed.x);
        double a1 = std::atan2(v.point.y - reg.seed.y, v.point.x - reg.seed.x);
        if (a1 < a0) {
            a1 += 2.0 * M_PI;
        }
        double a_diff = a1 - a0;
        int angle_n = static_cast<int>(std::ceil(a_diff / max_sample_beta));
        double angle_delta = a_diff / static_cast<double>(angle_n);
        double alpha = a0;
        for (int j = 0; j < angle_n - 1; ++j) {
            alpha += angle_delta;
            FPoint p = MakePoint(reg.seed.x + reg.radius * std::cos(alpha), reg.seed.y + reg.radius * std::sin(alpha));
            ret.vertices.push_back({-1, -4, p});
        }
        ret.vertices.push_back(v);
        ret.vertices.back().edge_flag = -4;
    }
    return ret;
}
