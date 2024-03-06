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
    const trivis::Trivis &vis
) {
    drawer.DrawMap();
    drawer.DrawPlane(kColorBlack);
    drawer.DrawBorders(kColorWhite);
    drawer.DrawHoles(kColorDimGray);
    drawer.DrawPolygons(vis.triangles(), 0.01, kColorYellow);
    for (int i = 0; i < vis.mesh().vertices.size(); ++i) {
        drawer.DrawText(std::to_string(i), vis.mesh().point(i), 0.03, kColorNavy);
    }
    for (int i = 0; i < vis.mesh().edges.size(); ++i) {
        drawer.DrawText(std::to_string(i), (vis.mesh().point(vis.mesh().edges[i].vertices[0]) + vis.mesh().point(vis.mesh().edges[i].vertices[1])) / 2.0, 0.03, kColorMaroon);
    }
    for (int i = 0; i < vis.mesh().triangles.size(); ++i) {
        const auto &tri = vis.mesh().triangles[i];
        drawer.DrawText(std::to_string(i), (vis.mesh().point(tri.vertices[0]) + vis.mesh().point(tri.vertices[1]) + vis.mesh().point(tri.vertices[2])) / 3.0, 0.03, kColorDeepPink);
    }
}

void drawing::FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const VisibilityRegionColors &colors,
    double w,
    double s
) {
    drawer.DrawArc(vis_reg.seed, 0.60, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawArc(vis_reg.seed, 0.40, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawArc(vis_reg.seed, 0.20, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawArc(vis_reg.seed, 0.16, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawArc(vis_reg.seed, 0.12, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawArc(vis_reg.seed, 0.08, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawArc(vis_reg.seed, 0.04, 0.0, 2.0 * M_PI, w * 0.02, colors.point);
    drawer.DrawPoint(vis_reg.seed, w * 0.02, colors.point);
    auto reg_with_arcs_samples = vis_reg;
    reg_with_arcs_samples.SampleArcEdges(M_PI / 180);
    auto polygon = reg_with_arcs_samples.ToPolygon();
    double vis_radius = vis_reg.radius.value_or(0.0);
    if (vis_radius > 0) {
        drawer.DrawArc(vis_reg.seed, vis_radius, 0.0, 2.0 * M_PI, w * 0.3, colors.base, 0.2);
    }
    drawer.DrawPolygon(polygon, colors.base, 0.5);
    if (vis_reg.vertices.size() <= 1) {
        drawer.DrawArc(vis_reg.seed, vis_radius, 0.0, 2.0 * M_PI, w * 0.3, colors.edge_arc);
        return;
    }

    for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
        const auto &vi_prev = vis_reg.vertices[i_prev];
        const auto &vi = vis_reg.vertices[i];
        if (vi.edge_flag > -1) {
            drawer.DrawLine(vi_prev.point, vi.point, w * 0.2, colors.edge_obstacle);
        } else if (vi.edge_flag == -1) {
            drawer.DrawLine(vi_prev.point, vi.point, w * 0.2, colors.edge_free);
        } else if (vi.edge_flag == -2) {
            drawer.DrawLine(vi_prev.point, vi.point, w * 0.2, colors.base, 0.2);
            auto a = vi.point - vis_reg.seed;
            auto b = vi_prev.point - vis_reg.seed;
            drawer.DrawArc(vis_reg.seed, vis_radius, std::atan2(a.x, a.y) - M_PI_2, std::atan2(b.x, b.y) - M_PI_2, w * 0.2, colors.edge_arc);
        } else if (vi.edge_flag == -3) {
            drawer.DrawLine(vi_prev.point, vi.point, w * 0.2, colors.edge_other);
        } else {
            drawer.DrawLine(vi_prev.point, vi.point, w * 0.2, colors.edge_other);
        }
    }
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
        drawer.DrawText(text, vi.point + u * norm * 0.5, s * 0.06, kColorBlack);
    }
    for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
        const auto &vi_prev = vis_reg.vertices[i_prev];
        const auto &vi = vis_reg.vertices[i];
        if (vi.edge_flag > 0) {
            drawer.DrawPoint(vi.point, w * 0.14, colors.edge_obstacle);
        } else if (vi.edge_flag == -1) {
            drawer.DrawPoint(vi.point, w * 0.14, colors.edge_free);
        } else if (vi.edge_flag == -2) {
            drawer.DrawPoint(vi.point, w * 0.14, colors.edge_arc);
        } else if (vi.edge_flag == -3) {
            drawer.DrawPoint(vi.point, w * 0.14, colors.edge_other);
        } else {
            drawer.DrawPoint(vi.point, w * 0.14, colors.edge_other);
        }
        drawer.DrawPoint(vi.point, w * 0.10, vi.vertex_flag == -1 ? colors.vertex_intersection : colors.vertex_node);
    }
    for (int i_prev = static_cast<int>(vis_reg.vertices.size()) - 1, i = 0; i < vis_reg.vertices.size(); i_prev = i++) {
        const auto &vi_prev = vis_reg.vertices[i_prev];
        const auto &vi = vis_reg.vertices[i];
        std::string text = std::to_string(i);
        if (vi.vertex_flag >= 0) {
            text += "-" + std::to_string(vi.vertex_flag);
        }
        drawer.DrawText(text, vi.point, s * 0.06, kColorBlack);
    }
}

void drawing::FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const RGB &c,
    double w,
    double s
) {
    RGB c_dark = {static_cast<int>(c.r / 1.5), static_cast<int>(c.g / 1.5), static_cast<int>(c.b / 1.5)};
    RGB c_darker = {static_cast<int>(c.r / 3.0), static_cast<int>(c.g / 3.0), static_cast<int>(c.b / 3.0)};
    RGB c_darkest = {static_cast<int>(c.r / 4.5), static_cast<int>(c.g / 4.5), static_cast<int>(c.b / 4.5)};
    VisibilityRegionColors colors;
    colors.base = c;
    colors.point = c_darker;
    colors.edge_other = c_darkest;
    colors.edge_obstacle = c_darker;
    colors.edge_free = c_dark;
    colors.edge_arc = c;
    colors.vertex_node = kColorWhite;
    colors.vertex_intersection = c_dark;
    FancyDrawRadialVisibilityRegion(drawer, vis_reg, colors, w, s);
}