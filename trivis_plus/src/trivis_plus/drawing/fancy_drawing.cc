/**
 * File:   fancy_drawing.cc
 *
 * Date:   08.11.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/drawing/fancy_drawing.h"

using namespace trivis_plus;
using namespace trivis_plus::drawing;

using namespace trivis;
using namespace trivis::geom;

void drawing::FancyDrawMap(
    const MapDrawer &drawer,
    const trivis::Trivis &vis,
    const MapColors &colors,
    double line_width,
    bool print_labels,
    double label_scale
) {
    drawer.DrawMap();
    drawer.DrawPlane(colors.plane);
    drawer.DrawBorders(colors.borders);
    drawer.DrawHoles(colors.holes);
    drawer.DrawPolygons(vis.triangles(), line_width * 0.05, colors.triangles);
    if (print_labels) {
        for (int i = 0; i < vis.mesh().vertices.size(); ++i) {
            drawer.DrawText(std::to_string(i), vis.mesh().point(i), label_scale * 0.05, kColorNavy);
        }
        for (int i = 0; i < vis.mesh().edges.size(); ++i) {
            auto midpoint = (vis.mesh().point(vis.mesh().edges[i].vertices[0]) + vis.mesh().point(vis.mesh().edges[i].vertices[1])) / 2.0;
            drawer.DrawText(std::to_string(i), midpoint, label_scale * 0.05, kColorMaroon);
        }
        for (int i = 0; i < vis.mesh().triangles.size(); ++i) {
            const auto &tri = vis.mesh().triangles[i];
            auto midpoint = (vis.mesh().point(tri.vertices[0]) + vis.mesh().point(tri.vertices[1]) + vis.mesh().point(tri.vertices[2])) / 3.0;
            drawer.DrawText(std::to_string(i), midpoint, label_scale * 0.05, kColorDeepPink);
        }
    }
}

void drawing::FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const VisibilityRegionColors &colors,
    double line_width,
    bool print_labels,
    double label_scale
) {
    auto reg_with_arcs_samples = vis_reg;
    reg_with_arcs_samples.SampleArcEdges(M_PI / 180);
    auto polygon = reg_with_arcs_samples.ToPolygon();
    double vis_radius = vis_reg.radius.value_or(0.0);
    if (vis_radius > 0) {
        drawer.DrawArc(vis_reg.seed, vis_radius, 0.0, 2.0 * M_PI, line_width * 0.2, colors.edge_arc, colors.edge_opacity * 0.2);
    }
    drawer.DrawPolygon(polygon, colors.base, colors.base_opacity);
    if (vis_reg.vertices.size() <= 1) {
        drawer.DrawArc(vis_reg.seed, vis_radius, 0.0, 2.0 * M_PI, line_width * 0.2, colors.edge_arc, colors.edge_opacity);
        return;
    }
    for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
        const auto &vi_prev = vis_reg.vertices[i_prev];
        const auto &vi = vis_reg.vertices[i];
        if (vi.edge_flag > -1) {
            drawer.DrawLine(vi_prev.point, vi.point, line_width * 0.2, colors.edge_obstacle, colors.edge_opacity);
        } else if (vi.edge_flag == -1) {
            drawer.DrawLine(vi_prev.point, vi.point, line_width * 0.2, colors.edge_free, colors.edge_opacity);
        } else if (vi.edge_flag == -2) {
            drawer.DrawLine(vi_prev.point, vi.point, line_width * 0.2, colors.base, colors.edge_opacity * 0.2);
            auto a = vi.point - vis_reg.seed;
            auto b = vi_prev.point - vis_reg.seed;
            drawer.DrawArc(vis_reg.seed, vis_radius, std::atan2(a.x, a.y) - M_PI_2, std::atan2(b.x, b.y) - M_PI_2, line_width * 0.2, colors.edge_arc, colors.edge_opacity);
        } else if (vi.edge_flag == -3) {
            drawer.DrawLine(vi_prev.point, vi.point, line_width * 0.2, colors.edge_sample, colors.edge_opacity);
        } else {
            drawer.DrawLine(vi_prev.point, vi.point, line_width * 0.2, colors.edge_unknown, colors.edge_opacity);
        }
    }
    drawer.DrawPoint(vis_reg.seed, line_width * 0.3, colors.point, colors.point_opacity);
    if (print_labels) {
        for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
            const auto &vi_prev = vis_reg.vertices[i_prev];
            const auto &vi = vis_reg.vertices[i];
            auto u = vi_prev.point - vi.point;
            double norm = u.Norm();
            u = u / u.Norm();
            std::string text = std::to_string(i);
            if (vi.edge_flag >= 0) {
                text += "-" + std::to_string(vi.edge_flag);
            }
            drawer.DrawText(text, vi.point + u * norm * 0.5, label_scale * 0.05, kColorBlack);
        }
    }
    for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
        const auto &vi_prev = vis_reg.vertices[i_prev];
        const auto &vi = vis_reg.vertices[i];
        if (vi.edge_flag > 0) {
            drawer.DrawPoint(vi.point, line_width * 0.15, colors.edge_obstacle, colors.edge_opacity);
        } else if (vi.edge_flag == -1) {
            drawer.DrawPoint(vi.point, line_width * 0.15, colors.edge_free, colors.edge_opacity);
        } else if (vi.edge_flag == -2) {
            drawer.DrawPoint(vi.point, line_width * 0.15, colors.edge_arc, colors.edge_opacity);
        } else if (vi.edge_flag == -3) {
            drawer.DrawPoint(vi.point, line_width * 0.15, colors.edge_sample, colors.edge_opacity);
        } else {
            drawer.DrawPoint(vi.point, line_width * 0.15, colors.edge_unknown, colors.edge_opacity);
        }
        drawer.DrawPoint(vi.point, line_width * 0.10, vi.vertex_flag == -1 ? colors.vertex_intersection : colors.vertex_node, colors.edge_opacity);
    }
    if (print_labels) {
        for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
            const auto &vi_prev = vis_reg.vertices[i_prev];
            const auto &vi = vis_reg.vertices[i];
            std::string text = std::to_string(i);
            if (vi.vertex_flag >= 0) {
                text += "-" + std::to_string(vi.vertex_flag);
            }
            drawer.DrawText(text, vi.point, label_scale * 0.05, kColorBlack);
        }
    }
}

void drawing::FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const RGB &color,
    double line_width,
    bool print_labels,
    double label_scale
) {
    RGB c_dark = {static_cast<int>(color.r / 1.5), static_cast<int>(color.g / 1.5), static_cast<int>(color.b / 1.5)};
    RGB c_darker = {static_cast<int>(color.r / 3.0), static_cast<int>(color.g / 3.0), static_cast<int>(color.b / 3.0)};
    RGB c_darkest = {static_cast<int>(color.r / 4.5), static_cast<int>(color.g / 4.5), static_cast<int>(color.b / 4.5)};
    VisibilityRegionColors colors;
    colors.base = color;
    colors.point = c_darker;
    colors.edge_obstacle = c_darker;
    colors.edge_free = c_dark;
    colors.edge_arc = color;
    colors.edge_sample = c_darkest;
    colors.edge_unknown = c_darkest;
    colors.vertex_node = kColorWhite;
    colors.vertex_intersection = c_dark;
    FancyDrawRadialVisibilityRegion(drawer, vis_reg, colors, line_width, print_labels, label_scale);
}