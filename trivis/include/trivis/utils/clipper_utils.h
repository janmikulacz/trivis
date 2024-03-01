/**
 * File:   clipper_utils.h
 *
 * Date:   28.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_UTILS_CLIPPER_UTILS_H_
#define TRIVIS_UTILS_CLIPPER_UTILS_H_

#include "clipper2/clipper.h"
#include "trivis/geom/geom_types.h"

namespace trivis::utils {

static constexpr double kDefaultMultiplier = 1073741824.0; // 2^30
static constexpr double kDefaultNormalizer = 1.0;

inline double Number2Clipper(
    double x,
    double normalizer = kDefaultNormalizer,
    double multiplier = kDefaultMultiplier
) {
    return multiplier * (x / normalizer);
}

inline Clipper2Lib::PointD Geom2Clipper(
    const geom::FPoint &p,
    double normalizer = kDefaultNormalizer,
    const geom::FPoint &subtrahend = {0.0, 0.0},
    double multiplier = kDefaultMultiplier
) {
    return {multiplier * ((p.x - subtrahend.x) / normalizer), multiplier * ((p.y - subtrahend.y) / normalizer)};
}

inline geom::FPoint Clipper2Geom(
    const Clipper2Lib::PointD &p,
    double normalizer = kDefaultNormalizer,
    const geom::FPoint &subtrahend = {0.0, 0.0},
    double multiplier = kDefaultMultiplier
) {
    return {(p.x / multiplier) * normalizer + subtrahend.x, (p.y / multiplier) * normalizer + subtrahend.y};
}

Clipper2Lib::PathD Geom2Clipper(
    const geom::FPolygon &polygon,
    double normalizer = kDefaultNormalizer,
    const geom::FPoint &subtrahend = {0.0, 0.0},
    double multiplier = kDefaultMultiplier
);

inline double NumberFromClipper(
    double x,
    double normalizer = kDefaultNormalizer,
    double multiplier = kDefaultMultiplier
) {
    return normalizer * (x / multiplier);
}

geom::FPolygon Clipper2Geom(
    const Clipper2Lib::PathD &polygon,
    double normalizer = kDefaultNormalizer,
    const geom::FPoint &subtrahend = {0.0, 0.0},
    double multiplier = kDefaultMultiplier
);

Clipper2Lib::PathsD Geom2Clipper(
    const geom::FPolygons &polygons,
    double normalizer = kDefaultNormalizer,
    const geom::FPoint &subtrahend = {0.0, 0.0},
    double multiplier = kDefaultMultiplier
);

geom::FPolygons Clipper2Geom(
    const Clipper2Lib::PathsD &polygons,
    double normalizer = kDefaultNormalizer,
    const geom::FPoint &subtrahend = {0.0, 0.0},
    double multiplier = kDefaultMultiplier
);

}

#endif //TRIVIS_UTILS_CLIPPER_UTILS_H_
