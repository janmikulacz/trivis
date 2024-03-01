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

void FancyDrawMap(
    const MapDrawer &drawer,
    const trivis::Trivis &vis
);

struct VisibilityRegionColors {
    RGB point = kColorDeepSkyBlue;
    RGB base = kColorLightSkyBlue;
    RGB edge_obstacle = kColorRed;
    RGB edge_free = kColorLimeGreen;
    RGB edge_arc = kColorBlueViolet;
    RGB edge_other = kColorDeepPink;
    RGB vertex_node = kColorOrange;
    RGB vertex_intersection = kColorYellow;
};

void FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const VisibilityRegionColors &colors = VisibilityRegionColors{},
    double w = 1.0,
    double s = 1.0
);

void FancyDrawRadialVisibilityRegion(
    const MapDrawer &drawer,
    const trivis::RadialVisibilityRegion &vis_reg,
    const RGB &c,
    double w = 1.0,
    double s = 1.0
);

}

#endif //TRIVIS_PLUS_DRAWING_FANCY_DRAWING_H_
