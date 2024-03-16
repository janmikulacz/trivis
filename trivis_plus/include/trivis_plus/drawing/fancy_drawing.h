/**
 * File:   fancy_drawing.h
 *
 * Date:   08.11.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DRAWING_FANCY_DRAWING_H_
#define TRIVIS_PLUS_DRAWING_FANCY_DRAWING_H_

#include "trivis_plus/drawing/drawing.h"

#include "trivis/trivis.h"

namespace trivis_plus::drawing {

struct MapColors {
    RGB plane = kColorBlack;
    RGB borders = kColorWhite;
    RGB holes = kColorDimGray;
    RGB triangles = kColorYellow;
};

void FancyDrawMap(
    const MapDrawer &drawer,
    const trivis::Trivis &vis,
    const MapColors &colors = MapColors{},
    double line_width = 1.0,
    bool print_labels = false,
    double label_scale = 1.0
);

struct VisibilityRegionColors {
    double base_opacity = 0.5;
    double edge_opacity = 1.0;
    double point_opacity = 1.0;
    RGB point = kColorBlueViolet;
    RGB base = kColorLimeGreen;
    RGB edge_obstacle = kColorRed;
    RGB edge_free = kColorLimeGreen;
    RGB edge_arc = kColorBlueViolet;
    RGB edge_sample = kColorBlue;
    RGB edge_unknown = kColorDeepPink;
    RGB vertex_node = kColorOrange;
    RGB vertex_intersection = kColorYellow;
};

void FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const VisibilityRegionColors &colors = VisibilityRegionColors{},
    double line_width = 1.0,
    bool print_labels = false,
    double label_scale = 1.0
);

void FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const RGB &color,
    double line_width = 1.0,
    bool print_labels = false,
    double label_scale = 1.0
);

}

#endif //TRIVIS_PLUS_DRAWING_FANCY_DRAWING_H_
