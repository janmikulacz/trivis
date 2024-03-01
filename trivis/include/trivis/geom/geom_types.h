/**
 * File:   geom_types.h
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_GEOM_GEOM_TYPES_H_
#define TRIVIS_GEOM_GEOM_TYPES_H_

#include "trivis/geom/generic_geom_types.h"

namespace trivis::geom {

struct FLimits {
    double x_min = 0.0;
    double x_max = 1.0;
    double y_min = 0.0;
    double y_max = 1.0;
};

using FPoint = Point<double>;
using FPosition = Position<double>;
using FPositions = Positions<double>;
using FPoints = Points<double>;
using FPolygon = Polygon<double>;
using FPolygons = Polygons<double>;

}

#endif //TRIVIS_GEOM_GEOM_TYPES_H_
