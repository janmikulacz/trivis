/**
 * File:   drawing.cc
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/drawing/drawing.h"

using namespace trivis_plus;
using namespace trivis_plus::drawing;

using namespace trivis;
using namespace trivis::geom;

MapDrawer drawing::MakeMapDrawer(
    const FPolygons &borders,
    const FPolygons &holes,
    double res,
    double relative_frame_width
) {
    double x_min, x_max, y_min, y_max;
    ComputeLimits(borders, x_min, x_max, y_min, y_max);
    double off = relative_frame_width * std::max(x_max - x_min, y_max - y_min);
    return MapDrawer(borders, holes, x_max - x_min + 2.0 * off, y_max - y_min + 2.0 * off, 1.0 / res, -x_min + off, -y_min + off);
}
